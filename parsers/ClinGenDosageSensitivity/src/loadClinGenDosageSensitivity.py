import os
import enum
import gzip

from Common.extractor import Extractor
from Common.loader_interface import SourceDataLoader
from Common.node_types import PRIMARY_KNOWLEDGE_SOURCE
from Common.prefixes import HGNC  # only an example, use existing curie prefixes or add your own to the prefixes file
from Common.utils import GetData
from datetime import date


# if parsing a tsv or csv type file with columns, use a enum to represent each field
class ClinGenDosageSensitivityGeneCOLS(enum.IntEnum):
    # HI: HaploInsufficiency ; TS : TriploSensitivity
    GENE = 1
    HI_DISEASE_GENE = -2
    TS_DISEASE_GENE = -1
    HI_SCORE_GENE = 4
    HI_DESCRIPTION_GENE = 5
    TS_SCORE_GENE = 12
    TS_DESCRIPTION_GENE = 13

class ClinGenDosageSensitivityRegionCOLS(enum.IntEnum):
    # HI: HaploInsufficiency ; TS : TriploSensitivity
    REGION = 0
    HI_DISEASE_REGION = -2
    TS_DISEASE_REGION = -1
    HI_SCORE_REGION = 4
    HI_DESCRIPTION_REGION = 5
    TS_SCORE_REGION = 12
    TS_DESCRIPTION_REGION = 13



##############
# Class: XXXX source loader
#
# Desc: Class that loads/parses the XXXX data.
##############
class ClinGenDosageSensitivityLoader(SourceDataLoader):
    source_id: str = 'ClinGenDosageSensitivity'
    # this should be a valid infores curie from the biolink infores catalog
    provenance_id: str = 'infores:clingen'  # Need to figure out this # Can only be filled from one of the values from https://github.com/biolink/biolink-model/blob/master/infores_catalog.yaml
    # increment parsing_version whenever changes are made to the parser that would result in changes to parsing output
    parsing_version: str = 'v1.0'

    def __init__(self, test_mode: bool = False, source_data_dir: str = None):
        """
        :param test_mode - sets the run into test mode
        :param source_data_dir - the specific storage directory to save files in
        """
        super().__init__(test_mode=test_mode, source_data_dir=source_data_dir)

        self.cligen_dosage_sensitivity_url = 'ftp://ftp.clinicalgenome.org/'
<<<<<<< HEAD
        self.cligen_dosage_sensitivity_gene_file = 'ClinGen_gene_curation_list_GRCh38.tsv'
        self.clingen_dosage_sensitivity_region_file = "ClinGen_region_curation_list_GRCh38.tsv"
        self.data_files = [self.cligen_dosage_sensitivity_gene_file, self.clingen_dosage_sensitivity_region_file]
=======
        self.cligen_dosage_sensitivity_file = 'ClinGen_gene_curation_list_GRCh38.tsv'
        self.data_files = [self.cligen_dosage_sensitivity_file]
>>>>>>> 6ac21f3 (Added ClinGenDosageSensitivity Parser)

    def get_latest_source_version(self) -> str:
        # if possible go to the source and retrieve a string that is the latest version of the source data
        latest_version = date.today().strftime("%Y%m%d")
        return latest_version

    def get_data(self) -> bool:
        # get_data is responsible for fetching the files in self.data_files and saving them to self.data_path
<<<<<<< HEAD
        data_puller = GetData()
        for source in self.data_files:
            source_data_url = f'{self.cligen_dosage_sensitivity_url}{source}'
            data_puller.pull_via_http(source_data_url, self.data_path)
=======
        source_data_url = f'{self.cligen_dosage_sensitivity_url}{self.cligen_dosage_sensitivity_file}'
        data_puller = GetData()
        data_puller.pull_via_http(source_data_url, self.data_path)
>>>>>>> 6ac21f3 (Added ClinGenDosageSensitivity Parser)
        return True

    def parse_data(self) -> dict:
        """
        Parses the data file for graph nodes/edges

        :return: ret_val: load_metadata
        """
        # This is a made up example of how one might extract nodes and edges from a tsv file
        # In this case it's taking the subject ID from column 1 and the object ID from column 3,
        # prepending them with a curie prefix. The predicate comes from column 2. The value in column 4
        # is set as a property on the edge.
        extractor = Extractor(file_writer=self.output_file_writer)
<<<<<<< HEAD
        dosage_sensitivity_gene_file: str = os.path.join(self.data_path, self.cligen_dosage_sensitivity_gene_file)
        with open(dosage_sensitivity_gene_file, 'rt') as fp:
            extractor.csv_extract(fp,
                                  lambda line: f'NCBIGene:{line[ClinGenDosageSensitivityGeneCOLS.GENE.value]}',
                                  # subject id
                                  lambda line: 'MONDO:0700096' if line[ClinGenDosageSensitivityGeneCOLS.HI_DISEASE_GENE.value] == '' else line[ClinGenDosageSensitivityGeneCOLS.HI_DISEASE_GENE.value],  # object id
=======
        dosage_sensitivity_file: str = os.path.join(self.data_path, self.cligen_dosage_sensitivity_file)
        with open(dosage_sensitivity_file, 'rt') as fp:
            extractor.csv_extract(fp,
                                  lambda line: f'NCBIGene:{line[ClinGenDosageSensitivityCOLS.GENE.value]}',
                                  # subject id
                                  lambda line: f'{line[ClinGenDosageSensitivityCOLS.HI_DISEASE.value]}',  # object id
>>>>>>> 6ac21f3 (Added ClinGenDosageSensitivity Parser)
                                  lambda line: 'gene associated with condition',  # predicate extractor
                                  lambda line: {},  # subject properties
                                  lambda line: {},  # object properties
                                  lambda line: dict({PRIMARY_KNOWLEDGE_SOURCE: self.provenance_id,
                                                     'Haploinsufficiency_Description': line[
<<<<<<< HEAD
                                                         ClinGenDosageSensitivityGeneCOLS.HI_DESCRIPTION_GENE.value]},
                                                    **get_edge_properties(
                                                        line[ClinGenDosageSensitivityGeneCOLS.HI_SCORE_GENE.value], "Haploinsufficiency", line[ClinGenDosageSensitivityGeneCOLS.HI_DISEASE_GENE.value])),
=======
                                                         ClinGenDosageSensitivityCOLS.HI_DESCRIPTION.value]},
                                                    **get_edge_properties(
                                                        line[ClinGenDosageSensitivityCOLS.HI_SCORE.value], "Haploinsufficiency")),
>>>>>>> 6ac21f3 (Added ClinGenDosageSensitivity Parser)

                                  # edge properties
                                  comment_character='#',
                                  delim='\t',
                                  has_header_row=True)

<<<<<<< HEAD
        with open(dosage_sensitivity_gene_file, 'rt') as fp:
            extractor.csv_extract(fp,
                                  lambda line: f'NCBIGene:{line[ClinGenDosageSensitivityGeneCOLS.GENE.value]}',
                                  # subject id
                                  lambda line: 'MONDO:0700096' if line[ClinGenDosageSensitivityGeneCOLS.TS_DISEASE_GENE.value] == '' else line[ClinGenDosageSensitivityGeneCOLS.TS_DISEASE_GENE.value],  # object id
=======
        with open(dosage_sensitivity_file, 'rt') as fp:
            extractor.csv_extract(fp,
                                  lambda line: f'NCBIGene:{line[ClinGenDosageSensitivityCOLS.GENE.value]}',
                                  # subject id
                                  lambda line: f'{line[ClinGenDosageSensitivityCOLS.TS_DISEASE.value]}',  # object id
>>>>>>> 6ac21f3 (Added ClinGenDosageSensitivity Parser)
                                  lambda line: 'gene associated with condition',  # predicate extractor
                                  lambda line: {},  # subject properties
                                  lambda line: {},  # object properties
                                  lambda line: dict({PRIMARY_KNOWLEDGE_SOURCE: self.provenance_id,
                                                     'Triplosensitivity_Description': line[
<<<<<<< HEAD
                                                         ClinGenDosageSensitivityGeneCOLS.TS_DESCRIPTION_GENE.value]},
                                                    **get_edge_properties(
                                                        line[ClinGenDosageSensitivityGeneCOLS.TS_SCORE_GENE.value],"Triplosensitivity", line[ClinGenDosageSensitivityGeneCOLS.TS_DISEASE_GENE])),
                                  # edge properties
                                  comment_character='#',
                                  delim='\t',
                                  has_header_row=True)
        dosage_sensitivity_region_file: str = os.path.join(self.data_path, self.clingen_dosage_sensitivity_region_file)
        with open(dosage_sensitivity_region_file, 'rt') as fp:
            extractor.csv_extract(fp,
                                  lambda line: f'ISCARegion:{line[ClinGenDosageSensitivityRegionCOLS.REGION.value]}',
                                  # subject id
                                  lambda line: 'MONDO:0700096' if line[ClinGenDosageSensitivityRegionCOLS.HI_DISEASE_REGION.value] == '' else line[ClinGenDosageSensitivityRegionCOLS.HI_DISEASE_REGION.value],  # object id,  # object id
                                  lambda line: 'region associated with condition',  # predicate extractor
                                  lambda line: {},  # subject properties
                                  lambda line: {},  # object properties
                                  lambda line: dict({PRIMARY_KNOWLEDGE_SOURCE: self.provenance_id,
                                                     'Haploinsufficiency_Description': line[
                                                         ClinGenDosageSensitivityRegionCOLS.HI_DESCRIPTION_REGION.value]},
                                                    **get_edge_properties(
                                                        line[ClinGenDosageSensitivityRegionCOLS.HI_SCORE_REGION.value], "Haploinsufficiency", line[ClinGenDosageSensitivityRegionCOLS.HI_DISEASE_REGION])),

                                  # edge properties
                                  comment_character='#',
                                  delim='\t',
                                  has_header_row=True)

        with open(dosage_sensitivity_region_file, 'rt') as fp:
            extractor.csv_extract(fp,
                                  lambda line: f'ISCARegion:{line[ClinGenDosageSensitivityRegionCOLS.REGION.value]}',
                                  # subject id
                                  lambda line: 'MONDO:0700096' if line[ClinGenDosageSensitivityRegionCOLS.TS_DISEASE_REGION.value] == '' else line[ClinGenDosageSensitivityRegionCOLS.TS_DISEASE_REGION.value],  # object id,  # object id
                                  lambda line: 'region associated with condition',  # predicate extractor
                                  lambda line: {},  # subject properties
                                  lambda line: {},  # object properties
                                  lambda line: dict({PRIMARY_KNOWLEDGE_SOURCE: self.provenance_id,
                                                     'Triplosensitivity_Description': line[
                                                         ClinGenDosageSensitivityRegionCOLS.TS_DESCRIPTION_REGION.value]},
                                                    **get_edge_properties(
                                                        line[ClinGenDosageSensitivityRegionCOLS.TS_SCORE_REGION.value],"Triplosensitivity", line[ClinGenDosageSensitivityRegionCOLS.TS_DISEASE_REGION])),
=======
                                                         ClinGenDosageSensitivityCOLS.TS_DESCRIPTION.value]},
                                                    **get_edge_properties(
                                                        line[ClinGenDosageSensitivityCOLS.TS_SCORE.value],"Triplosensitivity")),
>>>>>>> 6ac21f3 (Added ClinGenDosageSensitivity Parser)
                                  # edge properties
                                  comment_character='#',
                                  delim='\t',
                                  has_header_row=True)
        return extractor.load_metadata


<<<<<<< HEAD

# Created a common function to take the score value and return attributes,
# this may have problems with TS and HI score where the mode of inheritance is captured
# need to figure the proper naming conventions
def get_edge_properties(score_value, disease_type, mondo_id):
    if mondo_id != '':
        try:
            score = int(score_value)
        except ValueError:
            return {"Status": "Not yet evaluated"}
        if score == 0:
            return {'negated': False, f"{disease_type}_Score": score}
        elif score in (1, 2, 3):
            return {"negated": False, f"{disease_type}_Score": score}
        elif score == 30:
            return {"negated": False, f"{disease_type}_Score": score, "mode_of_inheritence": "Autosomal Recessive"}
        elif score == 40:
            return {"negated": True, f"{disease_type}_Score": score}
    else:
        return {"negated":True}

=======
# Created a common function to take the score value and return attributes,
# this may have problems with TS and HI score where the mode of inheritance is captured
# need to figure the proper naming conventions
def get_edge_properties(score_value, disease_type):
    try:
        score = int(score_value)
    except ValueError:
        return {"Status": "Not yet evaluated"}
    if score == 0:
        return {'negated': False, f"{disease_type}_Score": score}
    elif score in (1, 2, 3):
        return {"negated": False, f"{disease_type}_Score": score}
    elif score == 30:
        return {"negated": False, f"{disease_type}_Score": score, "mode_of_inheritence": "Autosomal Recessive"}
    elif score == 40:
        return {"negated": True, f"{disease_type}_Score": score}
>>>>>>> 6ac21f3 (Added ClinGenDosageSensitivity Parser)
