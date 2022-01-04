#!/usr/bin/env python
# coding: utf-8

# 
# <h1 style="color:#fc6b03">Mapping metadata and ontology</h1>
# 
# ### Data and metadata
# - Data from <span style="color:orange">Young modulus</span> experiment of concrete material in .dat file (table format)
# - Metadata was extracted from the header of the file (listed all metadata in a csv file)
# 
# ### Ontology
# - The <span style="color:orange">ontology BWMD_mid</span> as mid ontology
# - <span style="color:orange">Concrete_Ontology_rdf</span> as specific usecase in LeBeDigital
# 
# ### The workflow of the script
# 1. Installing and importing the necessary libraries: <span style="color:orange">owlready, rdflib and pandas</span>
# 2. Using owlready to load the 2 ontology files and pandas to load csv file
# 3. Using rdflib to create a graph and add triples by adding data from csv file and ontology entities from ontology files
# 4. Exporting the graph as .ttl file (turtle file)
# 

# In[1]:

from owlready2 import *
import pandas as pd
import urllib.parse
from pathlib import Path
import os
import datetime


# In[2]:


baseDir0 = Path(__file__).resolve().parents[0]
baseDir1 = Path(__file__).resolve().parents[1]
baseDir2 = Path(__file__).resolve().parents[2]
ontologyPath = os.path.join(baseDir2,'ConcreteOntology')
metadataPath = os.path.join(baseDir0,'E-modul-processed-data/emodul_metadata.csv')
graphPath = os.path.join(baseDir0,'E-modul-processed-data/EM_Graph.ttl')
#importedOntologiesPath = os.path.join(baseDir2,'GraphCreation/Ontologies')
processedDataPath = os.path.join(os.path.join(os.path.join(baseDir0,'E-modul-processed-data'),'processeddata'),'processed_')

# <h3 style="color:#1f5dbf">load concrete material ontology for Emodul experiment</h3>   

# In[3]:


onto_path.append(".")
My_world = World()

cco_ontology = My_world.get_ontology(ontologyPath + "/MergedAllCoreOntology.owl").load()
mseo_mid = My_world.get_ontology(ontologyPath + "/MSEO_mid.owl").load()
PeriodicTable_ontology = My_world.get_ontology(ontologyPath + "/PeriodicTable.owl").load()
WCTmidonto = My_world.get_ontology(ontologyPath + "/WCTmid.owl").load()
CSTonto = My_world.get_ontology(ontologyPath + "/ConcreteStressTestOntologie.owl").load()
lebedigital_concrete = My_world.get_ontology(os.path.join(ontologyPath,'EM.xml')).load()
ConcreteMSEO_ontology = My_world.get_ontology(ontologyPath + "/Concrete_Ontology_MSEO.owl").load()

BWMD = My_world.get_namespace("https://www.materials.fraunhofer.de/ontologies/BWMD_ontology/mid#")
WCTmid = My_world.get_namespace("https://mobi.com/ontologies/6/2021/WCTmid#")
#WCT = My_world.get_namespace("https://mobi.com/ontologies/6/202https://mobi.com/ontologies/6/2021/WoodCompressionTest#1/WoodCompressionTest#")
CCO = My_world.get_namespace("http://www.ontologyrepository.com/CommonCoreOntologies/")
CST = My_world.get_namespace("https://mobi.com/ontologies/7/2021/ConcreteStressTestOntologie#")
CON = My_world.get_namespace('http://w3id.org/concrete')
OBO = My_world.get_namespace("http://purl.obolibrary.org/obo/")
My_world.get_namespace('https://purl.matolab.org/mseo/mid/')

lebedigital_concrete.imported_ontologies.append(cco_ontology)
lebedigital_concrete.imported_ontologies.append(PeriodicTable_ontology)
lebedigital_concrete.imported_ontologies.append(WCTmidonto)
lebedigital_concrete.imported_ontologies.append(CSTonto)
lebedigital_concrete.imported_ontologies.append(mseo_mid)




# <h3 style="color:#1f5dbf">metadata from Emodul experiments</h3>  

# In[7]:


data = pd.read_csv(metadataPath)



# In[9]:


data['experiment name'] = ['E-modul_experiment_' + data['experiment raw name'][i].replace(' ','_').replace('.','_')
                           for i in data.index
                          ]
data['sample name 1'] = [data['sample name'][i].replace(' ','_').replace('.','_')
                           for i in data.index
                          ]

# <h5 style="color:#1f5dbf">splitting string from remark in control and the value of control force/stress</h5>  

# In[11]:


data['control'] = [
    data['remark'][i].split()[0] for i in data.index
]
data['control value'] = [
    float(data['remark'][i].split()[1].replace(',','.')) for i in data.index
]
data['control unit'] = [
    data['remark'][i].split()[2] for i in data.index
]


# <h5 style="color:#1f5dbf">add 1 column for example file path in data</h5>  

# In[13]:


data['file path'] = [
    'https://github.com/BAMresearch/ModelCalibration/tree/main/usecases/Concrete/Data/E-modul'
    + '/' + data['experiment raw name'][i]
    for i in data.index
]

data['processed data file path'] = [
    processedDataPath
    + data['experiment raw name'][i].replace(' ','_').replace('.','_')
    + '.csv'
    for i in data.index
]


# <h5 style="color:#1f5dbf">convert string to number in the columns</h5>  

# In[15]:


data['weight_number'] = [
    float(data['weight'][i].replace(',','.'))
    for i in data.index
]
data['diameter_number'] = [
    float(data['diameter'][i].replace(',','.'))
    for i in data.index
]
data['length_number'] = [
    float(data['length'][i].replace(',','.'))
    for i in data.index
]





# <h3 style="color:#1f5dbf">Generate a graph and add triples to this graph</h3>  

# In[17]:


from rdflib import URIRef, Graph, Literal, BNode
from rdflib.namespace import FOAF, RDF
g = My_world.as_rdflib_graph()


# <h5 style="color:#1f5dbf">adding instances into classes</h5>  

# In[18]:


for i in data.index:
    with lebedigital_concrete:
        # add experiments as instances in class DeterminationOfSecantModulusOfElasticity
        g.add(
            (
                URIRef(urllib.parse.unquote(ConcreteMSEO_ontology.DeterminationOfSecantModulusOfElasticity(data['experiment name'][i].replace(' ','_').replace('.','_')).iri)), 
                RDF.type, 
                URIRef(urllib.parse.unquote(ConcreteMSEO_ontology.DeterminationOfSecantModulusOfElasticity.iri))
            )
        )
        # add specimens as instances in class Specimen
        g.add(
            (
                URIRef(urllib.parse.unquote(lebedigital_concrete.Specimen('Specimen_' + data['sample name 1'][i]).iri)), 
                RDF.type, 
                URIRef(urllib.parse.unquote(lebedigital_concrete.Specimen.iri))
            )
        )
        # add testers (prüfers) as instances in class Agent
        g.add(
            (
                URIRef(urllib.parse.unquote(CCO.Agent('Agent_' + data['sample name 1'][i] +data['tester'][i].replace(' ','_')).iri)), 
                RDF.type, 
                URIRef(urllib.parse.unquote(CCO.Agent.iri))
            )
        )
        # add person/tester in class DesignativeName
        g.add(
            (
                URIRef(urllib.parse.unquote(CCO.DesignativeName('DesignativeName_' + data['sample name 1'][i] + data['tester'][i].replace(' ','_')).iri)), 
                RDF.type, 
                URIRef(urllib.parse.unquote(CCO.DesignativeName.iri))
            )
        )
        # add start date as instances in class Day
        g.add(
            (
                URIRef(urllib.parse.unquote(CCO.Day('Day_' + data['data collection timestamp'][i].replace(' ','_')).iri)), 
                RDF.type, 
                URIRef(urllib.parse.unquote(CCO.Day.iri))
            )
        )
        # add controls as instances in class ForceRate
        g.add(
            (
                URIRef(urllib.parse.unquote(ConcreteMSEO_ontology.ForceRate('ForceRate_' + data['sample name 1'][i] + data['control'][i]).iri)), 
                RDF.type, 
                URIRef(urllib.parse.unquote(ConcreteMSEO_ontology.ForceRate.iri))
            )
        )
        # add unit in MeasurementUnitOfForceRate
        g.add(
            (
                URIRef(urllib.parse.unquote(CON.MeasurementUnitOfForceRate('MeasurementUnitOfForceRate_' + data['sample name 1'][i]+ data['control unit'][i]).iri)), 
                RDF.type, 
                URIRef(urllib.parse.unquote(CON.MeasurementUnitOfForceRate.iri))
            )
        )
        # add length of the specimen in class Length
        g.add(
            (
                URIRef(urllib.parse.unquote(CCO.Length('Length_' + data['sample name 1'][i] + data['length'][i].replace(',','_')).iri)), 
                RDF.type, 
                URIRef(urllib.parse.unquote(CCO.Length.iri))
            )
        )
        # add Diameter of the specimen in class Diameter
        g.add(
            (
                URIRef(urllib.parse.unquote(CCO.Diameter('Diameter_' + data['sample name 1'][i] + data['diameter'][i].replace(',','_')).iri)), 
                RDF.type, 
                URIRef(urllib.parse.unquote(CCO.Diameter.iri))
            )
        )
        # add Weight of the specimen in class Mass
        g.add(
            (
                URIRef(urllib.parse.unquote(CCO.Mass('Mass_' + data['sample name 1'][i] + data['weight'][i].replace(',','_')).iri)), 
                RDF.type, 
                URIRef(urllib.parse.unquote(CCO.Mass.iri))
            )
        )
#        
        # add specimen region as instances in class MeasurementRegion
        g.add(
            (
                URIRef(urllib.parse.unquote(lebedigital_concrete.Specimen('MeasurementRegion_' + data['sample name 1'][i]).iri)), 
                RDF.type, 
                URIRef(urllib.parse.unquote(lebedigital_concrete.MeasurementRegion.iri))
            )
        )
        # add weight, diameter, length, dataset path, tester name value, force rate value, day
        # in class InformationBearingEntity
        g.add(
            (
                URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + data['sample name 1'][i] + data['weight'][i].replace(',','_')).iri)), 
                RDF.type, 
                URIRef(urllib.parse.unquote(CCO.InformationBearingEntity.iri))
            )
        )
        g.add(
            (
                URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + data['sample name 1'][i] + data['diameter'][i].replace(',','_')).iri)), 
                RDF.type, 
                URIRef(urllib.parse.unquote(CCO.InformationBearingEntity.iri))
            )
        )
        g.add(
            (
                URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + data['sample name 1'][i] + data['length'][i].replace(',','_')).iri)), 
                RDF.type, 
                URIRef(urllib.parse.unquote(CCO.InformationBearingEntity.iri))
            )
        )
        g.add(
            (
                URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + data['sample name 1'][i] + 'specimen_dat').iri)), 
                RDF.type, 
                URIRef(urllib.parse.unquote(CCO.InformationBearingEntity.iri))
            )
        )
        g.add(
            (
                URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + data['tester'][i].replace(' ','_')).iri)), 
                RDF.type, 
                URIRef(urllib.parse.unquote(CCO.InformationBearingEntity.iri))
            )
        )
        g.add(
            (
                URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + data['sample name 1'][i] + '_' + data['control'][i] ).iri)), 
                RDF.type, 
                URIRef(urllib.parse.unquote(CCO.InformationBearingEntity.iri))
            )
        )
        g.add(
            (
                URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + data['data collection timestamp'][i].replace(' ','_')).iri)), 
                RDF.type, 
                URIRef(urllib.parse.unquote(CCO.InformationBearingEntity.iri))
            )
        )

        # add processed data into class bfo:BFO_0000015
        g.add(
            (
                URIRef(urllib.parse.unquote(OBO.BFO_0000015('BFO_0000015_' + 'processed_' + data['sample name 1'][i]).iri)), 
                RDF.type, 
                URIRef(urllib.parse.unquote(OBO.BFO_0000015.iri))
            )
        )
        # add processed data into class AnalysedDataSet
        g.add(
            (
                URIRef(urllib.parse.unquote(lebedigital_concrete.AnalysedDataSet('AnalysedDataSet_' + 'processed_' + data['sample name 1'][i]).iri)), 
                RDF.type, 
                URIRef(urllib.parse.unquote(lebedigital_concrete.AnalysedDataSet.iri))
            )
        )
        # add processed data into class InformationBearingEntity
        g.add(
            (
                URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + 'processed_' + data['sample name 1'][i]).iri)), 
                RDF.type, 
                URIRef(urllib.parse.unquote(CCO.InformationBearingEntity.iri))
            )
        )


# <h5 style="color:#1f5dbf">adding data type in process dataset</h5>  

# In[19]:


for i in data.index:
    with lebedigital_concrete:
        g.add(
                (
                    URIRef(urllib.parse.unquote(lebedigital_concrete.RawDataSet('RawDataSet_' + data['sample name 1'][i] + 'specimen_dat').iri)), 
                    RDF.type, 
                    URIRef(urllib.parse.unquote(lebedigital_concrete.RawDataSet.iri))
                )
            )

# -------------------------------------------------------------------------------------------------------

# <h5 style="color:#1f5dbf">adding data with object properties and data properties</h5>  

# In[20]:


for i in data.index:
    with lebedigital_concrete:
        # experiment has_output RawDataSet ('specimen_dat')
        g.add(
            (
                URIRef(urllib.parse.unquote(ConcreteMSEO_ontology.DeterminationOfSecantModulusOfElasticity(data['experiment name'][i].replace(' ','_'))
                                          .iri)), 
                URIRef(urllib.parse.unquote(CCO.has_output.iri)), 
                URIRef(urllib.parse.unquote(lebedigital_concrete.RawDataSet('RawDataSet_' + data['sample name 1'][i] + 'specimen_dat').iri))
            )
        )
        # experiment has_agent from class person
        g.add(
            (
                URIRef(urllib.parse.unquote(ConcreteMSEO_ontology.DeterminationOfSecantModulusOfElasticity(data['experiment name'][i].replace(' ','_').replace('.','_'))
                                          .iri)), 
                URIRef(urllib.parse.unquote(CCO.has_agent.iri)), 
                URIRef(urllib.parse.unquote(CCO.Agent('Agent_' + data['sample name 1'][i] +data['tester'][i].replace(' ','_')).iri))
            )
        )
        # specimen is_input_of experiment
        g.add(
            (
                URIRef(urllib.parse.unquote(lebedigital_concrete.Specimen('Specimen_' + data['sample name 1'][i]).iri)), 
                URIRef(urllib.parse.unquote(CCO.is_input_of.iri)), 
                URIRef(urllib.parse.unquote(ConcreteMSEO_ontology.DeterminationOfSecantModulusOfElasticity(data['experiment name'][i].replace(' ','_').replace('.','_'))
                                          .iri))
            )
        )
        # ForceRate is_input_of experiment
        g.add(
            (
                URIRef(urllib.parse.unquote(ConcreteMSEO_ontology.ForceRate('ForceRate_' + data['sample name 1'][i] + data['control'][i]).iri)), 
                URIRef(urllib.parse.unquote(CCO.is_input_of.iri)), 
                URIRef(urllib.parse.unquote(ConcreteMSEO_ontology.DeterminationOfSecantModulusOfElasticity(data['experiment name'][i].replace(' ','_').replace('.','_'))
                                          .iri))
            )
        )
        # Experiment occures_on Day
        g.add(
            (
                URIRef(urllib.parse.unquote(ConcreteMSEO_ontology.DeterminationOfSecantModulusOfElasticity(data['experiment name'][i].replace(' ','_').replace('.','_'))
                                          .iri)), 
                URIRef(urllib.parse.unquote(CCO.occures_on.iri)), 
                URIRef(urllib.parse.unquote(CCO.Day('Day_' + data['data collection timestamp'][i].replace(' ','_')).iri))
            )
        )
        # specimen obo:BFO_0000051 MeasurementRegion
        g.add(
            (
                URIRef(urllib.parse.unquote(lebedigital_concrete.Specimen('Specimen_' + data['sample name 1'][i]).iri)), 
                URIRef(urllib.parse.unquote(OBO.BFO_0000051.iri)), 
                URIRef(urllib.parse.unquote(lebedigital_concrete.Specimen('MeasurementRegion_' + data['sample name 1'][i]).iri))
            )
        )
        # MeasurementRegion obo:BFO_0000086 Diameter, Length, Mass
        g.add(
            (
                URIRef(urllib.parse.unquote(lebedigital_concrete.Specimen('MeasurementRegion_' + data['sample name 1'][i]).iri)), 
                URIRef(urllib.parse.unquote(OBO.RO_0000086.iri)), 
                URIRef(urllib.parse.unquote(CCO.Diameter('Diameter_' + data['sample name 1'][i] + data['diameter'][i].replace(',','_')).iri))
            )
        )
        g.add(
            (
                URIRef(urllib.parse.unquote(lebedigital_concrete.Specimen('MeasurementRegion_' + data['sample name 1'][i]).iri)), 
                URIRef(urllib.parse.unquote(OBO.RO_0000086.iri)), 
                URIRef(urllib.parse.unquote(CCO.Length('Length_' + data['sample name 1'][i] + data['length'][i].replace(',','_')).iri))
            )
        )
        g.add(
            (
                URIRef(urllib.parse.unquote(lebedigital_concrete.Specimen('MeasurementRegion_' + data['sample name 1'][i]).iri)), 
                URIRef(urllib.parse.unquote(OBO.RO_0000086.iri)), 
                URIRef(urllib.parse.unquote(CCO.Mass('Mass_' + data['sample name 1'][i] + data['weight'][i].replace(',','_')).iri))
            )
        )
        # Diameter, Length, Mass, ForceRate, RawDataSet, DesignativeName obo:RO_0010001 InformationBearingEntity
        # Day obo:RO_0010001 InformationBearingEntity
        g.add(
            (
                URIRef(urllib.parse.unquote(CCO.Diameter('Diameter_' + data['sample name 1'][i] + data['diameter'][i].replace(',','_')).iri)), 
                URIRef(urllib.parse.unquote(OBO.RO_0010001.iri)), 
                URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + data['sample name 1'][i] + data['diameter'][i].replace(',','_')).iri))
            )
        )
        g.add(
            (
                URIRef(urllib.parse.unquote(CCO.Length('Length_' + data['sample name 1'][i] + data['length'][i].replace(',','_')).iri)), 
                URIRef(urllib.parse.unquote(OBO.RO_0010001.iri)), 
                URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + data['sample name 1'][i] + data['length'][i].replace(',','_')).iri))
            )
        )
        g.add(
            (
                URIRef(urllib.parse.unquote(CCO.Mass('Mass_' + data['sample name 1'][i] + data['weight'][i].replace(',','_')).iri)), 
                URIRef(urllib.parse.unquote(OBO.RO_0010001.iri)), 
                URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + data['sample name 1'][i] + data['weight'][i].replace(',','_')).iri))
            )
        )
        g.add(
            (
                URIRef(urllib.parse.unquote(ConcreteMSEO_ontology.ForceRate('ForceRate_' + data['sample name 1'][i] + data['control'][i]).iri)), 
                URIRef(urllib.parse.unquote(OBO.RO_0010001.iri)), 
                URIRef(urllib.parse.unquote(CCO.InformationBearingEntity(data['control'][i]).iri))
            )
        )
        g.add(
            (
                URIRef(urllib.parse.unquote(lebedigital_concrete.RawDataSet('RawDataSet_' + data['sample name 1'][i] + 'specimen_dat').iri)), 
                URIRef(urllib.parse.unquote(OBO.RO_0010001.iri)), 
                URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + data['sample name 1'][i] + 'specimen_dat').iri))
            )
        )
        g.add(
            (
                URIRef(urllib.parse.unquote(CCO.DesignativeName('DesignativeName_' + data['sample name 1'][i] + data['tester'][i].replace(' ','_')).iri)), 
                URIRef(urllib.parse.unquote(OBO.RO_0010001.iri)),
                URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + data['tester'][i].replace(' ','_')).iri))
            )
        )

        g.add(
            (
                URIRef(urllib.parse.unquote(CCO.Day('Day_' + data['data collection timestamp'][i].replace(' ','_')).iri)), 
                URIRef(urllib.parse.unquote(OBO.RO_0010001.iri)),
                URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + data['data collection timestamp'][i].replace(' ','_')).iri))
            )
        )
        # Agent designated_by DesignativeName
        g.add(
            (
                URIRef(urllib.parse.unquote(CCO.Agent('Agent_' + data['sample name 1'][i] +data['tester'][i].replace(' ','_')).iri)), 
                URIRef(urllib.parse.unquote(CCO.designated_by.iri)),
                URIRef(urllib.parse.unquote(CCO.DesignativeName('DesignativeName_' + data['sample name 1'][i] + data['tester'][i].replace(' ','_')).iri))
            )
        )
        # InformationBearingEntity of ForceRate uses_measurement_unit MeasurementUnitOfForceRate
        g.add(
            (
                URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + data['sample name 1'][i] + '_' + data['control'][i] ).iri)), 
                URIRef(urllib.parse.unquote(CCO.uses_measurement_unit.iri)), 
                URIRef(urllib.parse.unquote(CON.MeasurementUnitOfForceRate('MeasurementUnitOfForceRate_' + data['sample name 1'][i]+ data['control unit'][i]).iri))
            )
        )
        # RawDataSet cco:is_input_of bfo:BFO_0000015
        g.add(
                (
                    URIRef(urllib.parse.unquote(lebedigital_concrete.RawDataSet('RawDataSet_' + data['sample name 1'][i] + 'specimen_dat').iri)), 
                    URIRef(urllib.parse.unquote(CCO.is_input_of.iri)), 
                    URIRef(urllib.parse.unquote(OBO.BFO_0000015('BFO_0000015_' + 'processed_' + data['sample name 1'][i]).iri))
                )
            )

        # bfo:BFO_0000015 cco:has_output mseo:AnalysedDataSet
        g.add(
                (
                    URIRef(urllib.parse.unquote(OBO.BFO_0000015('BFO_0000015_' + 'processed_' + data['sample name 1'][i]).iri)), 
                    URIRef(urllib.parse.unquote(CCO.has_output.iri)), 
                    URIRef(urllib.parse.unquote(lebedigital_concrete.AnalysedDataSet('AnalysedDataSet_' + 'processed_' + data['sample name 1'][i]).iri))
                )
            )
        # mseo:AnalysedDataSet obo:RO_0010001 cco:InformationBearingEntity
        g.add(
                (
                    URIRef(urllib.parse.unquote(lebedigital_concrete.AnalysedDataSet('AnalysedDataSet_' + 'processed_' + data['sample name 1'][i]).iri)), 
                    URIRef(urllib.parse.unquote(OBO.RO_0010001.iri)), 
                    URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + 'processed_' + data['sample name 1'][i]).iri))
                )
            )

#-----------------------------------------------------------------------------------------------------

# <h5 style="color:#1f5dbf">adding data with data properties</h5>  

# In[21]:


# add weight, diameter, length, dataset path, tester name value, force rate value
for i in data.index:
    with lebedigital_concrete:
        g.add(
            (
                URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + data['sample name 1'][i] + data['weight'][i].replace(',','_')).iri)), 
                URIRef(urllib.parse.unquote(CCO.has_decimal_value.iri)),
                Literal(data['weight_number'][i])
            )
        )
        g.add(
            (
                URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + data['sample name 1'][i] + data['diameter'][i].replace(',','_')).iri)), 
                URIRef(urllib.parse.unquote(CCO.has_decimal_value.iri)),
                Literal(data['diameter_number'][i])
            )
        )
        g.add(
            (
                URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + data['sample name 1'][i] + data['length'][i].replace(',','_')).iri)), 
                URIRef(urllib.parse.unquote(CCO.has_decimal_value.iri)),
                Literal(data['length_number'][i])
            )
        )
        g.add(
            (
                URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + data['sample name 1'][i] + 'specimen_dat').iri)), 
                URIRef(urllib.parse.unquote(CCO.has_URI_value.iri)), 
                Literal(data['file path'][i])
            )
        )
        g.add(
            (
                URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + data['tester'][i].replace(' ','_')).iri)), 
                URIRef(urllib.parse.unquote(CCO.has_text_value.iri)), 
                Literal(data['tester'][i])
            )
        )
        g.add(
            (
                URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + data['sample name 1'][i] + '_' + data['control'][i]).iri)), 
                URIRef(urllib.parse.unquote(CCO.has_decimal_value.iri)), 
                Literal(data['control value'][i])
            )
        )
        g.add(
            (
                URIRef(urllib.parse.unquote(CON.MeasurementUnitOfForceRate('MeasurementUnitOfForceRate_' + data['sample name 1'][i]+ data['control unit'][i]).iri)), 
                URIRef(urllib.parse.unquote(CCO.has_text_value.iri)), 
                Literal(data['control unit'][i])
            )
        )
        # cco: InformationBearingEntity of processed dataset has filepath 
        g.add(
                (
                    URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + 'processed_' + data['sample name 1'][i]).iri)), 
                    URIRef(urllib.parse.unquote(CCO.has_URI_value.iri)), 
                    Literal(data['processed data file path'][i])
                )
            )

        # Day has_datetime_value
        g.add(
                (
                    URIRef(urllib.parse.unquote(CCO.InformationBearingEntity('InformationBearingEntity_' + data['data collection timestamp'][i].replace(' ','_')).iri)), 
                    URIRef(urllib.parse.unquote(CCO.has_datetime_value.iri)), 
                    Literal(datetime.datetime.strptime(data['data collection timestamp'][i], '%d.%m.%Y %H:%M:%S'))
                )
            )


g.serialize(destination=graphPath, format="turtle")

print('---------------------------------------------------------------------')
print('number of EModul experiments: ', data.shape[0])
print('number of classes in Emodul ontology: ', len(list(lebedigital_concrete.classes())))

q = """
    SELECT (count(*) as ?Triples)
    WHERE {
        {?s ?p ?o}
    }
"""


# In[31]:


results = g.query(q)
for result in results:
    print('number of triples', result)

