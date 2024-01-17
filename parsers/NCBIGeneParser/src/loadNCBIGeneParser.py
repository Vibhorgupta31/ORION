import os
import enum
import gzip

from Common.extractor import Extractor
from Common.kgxmodel import kgxnode, kgxedge
from Common.loader_interface import SourceDataLoader
from Common.node_types import PRIMARY_KNOWLEDGE_SOURCE
from Common.prefixes import HGNC  # only an example, use existing curie prefixes or add your own to the prefixes file
from Common.utils import GetData
import urllib.parse


##############
# Class: XXXX source loader
#
# Desc: Class that loads/parses the XXXX data.
##############

# http://useast.ensembl.org/info/website/upload/gff.html
class NCBIGeneParserCOLS(enum.IntEnum):
    SEQUENCE = 0
    TYPE = 2
    START_LOCATION = 3
    END_LOCTATION = 4
    STRAND = 6
    ATTRIBUTES = 8


class NCBIGeneParser(SourceDataLoader):
    source_id: str = 'NCBI'
    # this should be a valid infores curie from the biolink infores catalog
    provenance_id: str = 'infores:NCBI'
    # increment parsing_version whenever changes are made to the parser that would result in changes to parsing output
    parsing_version: str = '1.0'

    def __init__(self, test_mode: bool = False, source_data_dir: str = None):
        """
        :param test_mode - sets the run into test mode
        :param source_data_dir - the specific storage directory to save files in
        """
        super().__init__(test_mode=test_mode, source_data_dir=source_data_dir)

        self.ncbi_gene_data_url = 'https://ftp.ncbi.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.40_GRCh38.p14/'
        self.ncbi_gene_data_file = 'GCF_000001405.40_GRCh38.p14_genomic.gff.gz'
        self.data_files = [self.ncbi_gene_data_file]

    def get_latest_source_version(self) -> str:
        # if possible go to the source and retrieve a string that is the latest version of the source data
        latest_version = 'v1.0'
        return latest_version

    def get_data(self) -> bool:
        # get_data is responsible for fetching the files in self.data_files and saving them to self.data_path
        source_data_url = f'{self.ncbi_gene_data_url}{self.ncbi_gene_data_file}'
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
        # prepending them with a curie prefix. The predicate comes from column 3. The value in column 4
        # is set as a property on the edge.
        extractor = Extractor(file_writer=self.output_file_writer)
        ncbi_gene_data_file : str = os.path.join(self.data_path, self.ncbi_gene_data_file)
        with gzip.open(ncbi_gene_data_file, 'rt') as fp:
            extractor.csv_extract(fp,
                                  lambda line: f'NCBIGene:{gff_attributes2dict(line[NCBIGeneParserCOLS.ATTRIBUTES.value])[0]}',  # subject id # get_ncbi_id(line[GENERICDATACOLS.ATTRIBUTES.value])
                                  lambda line: f'HGNC:{gff_attributes2dict(line[NCBIGeneParserCOLS.ATTRIBUTES.value])[1]}', # get_hgnc_id(line[GENERICDATACOLS.ATTRIBUTES.value])
                                  # object id
                                  lambda line: "similar_to",  # predicate extractor
                                  lambda line: dict({'SEQUENCE': line[NCBIGeneParserCOLS.SEQUENCE.value], 'START_POSITION':line[NCBIGeneParserCOLS.START_LOCATION.value],
                                                     'END_POSITION': line[NCBIGeneParserCOLS.END_LOCTATION.value],
                                                     'STRAND':line[NCBIGeneParserCOLS.STRAND.value]}),
                                  # subject properties
                                  lambda line: {},  # object properties
                                  lambda line: {},  # edge properties
                                  comment_character='#',
                                  delim='\t',
                                  has_header_row=True,
                                  filter_set=["gene"],
                                  filter_field=NCBIGeneParserCOLS.TYPE.value)
        return extractor.load_metadata


def process_gff_kv(kv):
    k, v = kv.split('=')

    if ',' in v:

        # list of values...

        return (k, [urllib.parse.unquote(x) for x in v.split(',')])

    else:

        return (k, urllib.parse.unquote(v))


def gff_attributes2dict(attr):
    data_dict = dict(process_gff_kv(kv) for kv in attr.split(';'))
    dbxref = data_dict["Dbxref"]
    dbxref_dict = {}
    for items in dbxref:
        key = items.split(":")[0]
        value = items.split(":")[1:]
        dbxref_dict[key] = value
    if "GeneID" in dbxref_dict.keys() and "HGNC" in dbxref_dict.keys():
        gene_id = dbxref_dict["GeneID"][0]
        hgnc_id = dbxref_dict["HGNC"][-1]
        return (gene_id, hgnc_id)
    else:
        pass