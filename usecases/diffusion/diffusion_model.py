import fenics as fe
import numpy as np
import yaml
import sys
from pathlib import Path
from matplotlib import pyplot as plt


class artificial_data:
    """
    Function to generate artifical data based on available analytical solutions
    
    f_diff_plane_sheet_const : Fick diffusion into a plane sheet with uniform initial distribution and equal surface concentrations
    """
    def __init__(self):
        self.diff_coeff = None
        self.L = None
        self.C_L = None
        self.C_R = None
        self.C_0 = None
        self.nt = None
        self.dt = None
        self.n = None
        self.sigma = None

    def load_parameter(self, filename):
        '''
        Load parameters
    
        Parameter:
        ----------
        filename : string
            YAML file with parameters
    
        Returns:
        --------
        parameter: list
        '''
        with open(filename, 'r') as parameter:
            para = yaml.load(parameter, Loader=yaml.FullLoader)
        self.diff_coef = para['diff_coef']
        self.L = para['L']
        self.C_L = para['C_L']
        self.C_R = para['C_R']
        self.C_0 = para['C_0']
        self.nt = para['nt']
        self.dt = para['dt']
        self.n = para['n']
        self.sigma = para['sigma']
        return parameter

    def f_diff_plane_sheet_const(self, diff_coef, L, C_L, C_R, C_0, nt, dt, n, sigma):
            '''
            Fick diffusion into a plane sheet with
                - Uniform initial distribution
                - Equal surface concentrations
    
            Parameter:
            ----------
            diff_coef : float
                Diffusion coefficiant
            L : float
                Thickness sheet
            C_L : float
                Concentration left side
            C_R : float
                Concentration right side
            C_0 : float
                initial concentration
            nt : integer
                time steps
            dt : float
                step site (time)
            n : integer
                Addition count (Fouruer)
            sigma : float
                Standard deviation (noise)
    
            Returns:
            --------
            mass: list
                Amount of substance that has diffused into the area
            '''
            # Amount of substance diffused into the area after infinite time
            M_inf = L * (0.5 * (C_L + C_R) - C_0)
            #
            M = []
            for t_step in range(nt+1):
                M_t = M_inf
                for n_step in range(n+1):
                    k = (2*n_step + 1)**2 * np.pi**2
                    M_t = M_t - 8*M_inf/k * np.exp(-diff_coef * k * t_step * dt / L**2)
                M.append(M_t)
            # Generate some noise
            mu = 0
            noise = np.random.normal(mu, sigma, nt + 1)
            M_noise = []
            for i in range (nt + 1):
                M_noise.append(M[i] + noise[i])
            for i in range(10):
                if M_noise[i]<0:
                    M_noise[i]=0
            return M_noise


def main():
    '''
    Fick diffusion into a plane sheet with
        - Uniform initial distribution
        - Equal surface concentrations

    Parameter:
    ----------
    diff_coef : float
        Diffusion coefficiant

    Returns:
    --------
    mass: list
        Amount of substance that has diffused into the area
    '''
    metadata_files = Path(__file__).parent.glob('*_meta.yaml')
    #
    for file in metadata_files:
        virt_experiment = artificial_data()
        para = virt_experiment.load_parameter(file)
        mass_time = virt_experiment.f_diff_plane_sheet_const(virt_experiment.diff_coef, virt_experiment.L, virt_experiment.C_L, virt_experiment.C_R, virt_experiment.C_0, virt_experiment.nt, virt_experiment.dt, virt_experiment.n, virt_experiment.sigma)
        np.savetxt(Path(str(file).replace('meta.yaml', 'data.csv')), mass_time, delimiter =", ", fmt ='% s')
        plt.plot(mass_time)
        plt.savefig(Path(str(file).replace('meta.yaml', 'data.png')))

if __name__ == "__main__":
    #
    main()
