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
class ClinGenDosageSensitivityCOLS(enum.IntEnum):
    # HI: HaploInsufficiency ; TS : TriploSensitivity
    GENE = 1
    HI_DISEASE = -2
    TS_DISEASE = -1
    HI_SCORE = 4
    HI_DESCRIPTION = 5
    TS_SCORE = 12
    TS_DESCRIPTION = 13


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
        self.cligen_dosage_sensitivity_file = 'ClinGen_gene_curation_list_GRCh38.tsv'
        self.data_files = [self.cligen_dosage_sensitivity_file]

    def get_latest_source_version(self) -> str:
        # if possible go to the source and retrieve a string that is the latest version of the source data
        latest_version = date.today().strftime("%Y%m%d")
        return latest_version

    def get_data(self) -> bool:
        # get_data is responsible for fetching the files in self.data_files and saving them to self.data_path
        source_data_url = f'{self.cligen_dosage_sensitivity_url}{self.cligen_dosage_sensitivity_file}'
        data_puller = GetData()
        data_puller.pull_via_http(source_data_url, self.data_path)
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
        dosage_sensitivity_file: str = os.path.join(self.data_path, self.cligen_dosage_sensitivity_file)
        with open(dosage_sensitivity_file, 'rt') as fp:
            extractor.csv_extract(fp,
                                  lambda line: f'NCBIGene:{line[ClinGenDosageSensitivityCOLS.GENE.value]}',
                                  # subject id
                                  lambda line: f'{line[ClinGenDosageSensitivityCOLS.HI_DISEASE.value]}',  # object id
                                  lambda line: 'gene associated with condition',  # predicate extractor
                                  lambda line: {},  # subject properties
                                  lambda line: {},  # object properties
                                  lambda line: dict({PRIMARY_KNOWLEDGE_SOURCE: self.provenance_id,
                                                     'Haploinsufficiency_Description': line[
                                                         ClinGenDosageSensitivityCOLS.HI_DESCRIPTION.value]},
                                                    **get_edge_properties(
                                                        line[ClinGenDosageSensitivityCOLS.HI_SCORE.value], "Haploinsufficiency")),

                                  # edge properties
                                  comment_character='#',
                                  delim='\t',
                                  has_header_row=True)

        with open(dosage_sensitivity_file, 'rt') as fp:
            extractor.csv_extract(fp,
                                  lambda line: f'NCBIGene:{line[ClinGenDosageSensitivityCOLS.GENE.value]}',
                                  # subject id
                                  lambda line: f'{line[ClinGenDosageSensitivityCOLS.TS_DISEASE.value]}',  # object id
                                  lambda line: 'gene associated with condition',  # predicate extractor
                                  lambda line: {},  # subject properties
                                  lambda line: {},  # object properties
                                  lambda line: dict({PRIMARY_KNOWLEDGE_SOURCE: self.provenance_id,
                                                     'Triplosensitivity_Description': line[
                                                         ClinGenDosageSensitivityCOLS.TS_DESCRIPTION.value]},
                                                    **get_edge_properties(
                                                        line[ClinGenDosageSensitivityCOLS.TS_SCORE.value],"Triplosensitivity")),
                                  # edge properties
                                  comment_character='#',
                                  delim='\t',
                                  has_header_row=True)
        return extractor.load_metadata


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
