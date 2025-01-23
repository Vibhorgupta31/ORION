import csv
from typing import Mapping, Optional
from Common.kgxmodel import kgxnode, kgxedge
from Common.kgx_file_writer import KGXFileWriter
from Common.biolink_constants import PRIMARY_KNOWLEDGE_SOURCE, AGGREGATOR_KNOWLEDGE_SOURCES


class Extractor:
    """
    This is a class so that it can be used to accumulate nodes and edges across multiple files or input streams
    Also so that it can provide a few different interfaces (csv, sql) and keep the guts of the callback code in one
    place.
    """
    def __init__(self, file_writer: KGXFileWriter = None):
        # You might think it would be good to include all the extractors at this level, but they are really file or query
        # level things. You might want to use the same extractor with two differently formatted files or two different
        # sql queries.

        self.node_ids = set()
        self.nodes = []
        self.edges = []

        self.load_metadata = { 'record_counter': 0, 'skipped_record_counter': 0, 'errors': []}
        self.errors = []

        self.file_writer = file_writer

    def generator_extract(self, tuple_extractor):
        for edge in tuple_extractor:
            self.process_tuple(**edge)
    
    # TODO: test and see if this can replace original csv_extract
    def csv_extract2(self, infile,
                    subject_extractor,
                    object_extractor=None,
                    predicate_extractor=None,
                    subject_property_extractor=None,
                    object_property_extractor=None,
                    edge_property_extractor=None,
                    filter_set=None,
                    filter_field=None,
                    comment_character="#",
                    delim='\t',
                    has_header_row=False,
                    exclude_unconnected_nodes=False):
        """Read a csv, perform callbacks to retrieve node and edge info per row.
        Assumes that all of the properties extractable for a node occur on the line with the node identifier"""
        def tuple_extractor():
            for line in infile:
                if comment_character is not None and line.startswith(comment_character):
                    continue

                if has_header_row and not skipped_header:
                    skipped_header = True
                    continue

                if filter_field is not None:
                    filter_field_value = line[:-1].split(delim)[filter_field]
                    if filter_field_value not in filter_set:
                        continue

                reader = csv.reader([line], delimiter=delim)
                split_row = list(reader)[0]
                yield {
                    'subject_id': subject_extractor(split_row),
                    'object_id': object_extractor(split_row),
                    'predicate': predicate_extractor(split_row),
                    'subjectprops': subject_property_extractor(split_row),
                    'objectprops': object_property_extractor(split_row),
                    'edgeprops': edge_property_extractor(split_row),
                }
        self.generator_extract(tuple_extractor)

    def csv_extract(self, infile,
                    subject_extractor,
                    object_extractor=None,
                    predicate_extractor=None,
                    subject_property_extractor=None,
                    object_property_extractor=None,
                    edge_property_extractor=None,
                    filter_set=None,
                    filter_field=None,
                    comment_character="#",
                    delim='\t',
                    has_header_row=False,
                    exclude_unconnected_nodes=False):
        """Read a csv, perform callbacks to retrieve node and edge info per row.
        Assumes that all of the properties extractable for a node occur on the line with the node identifier"""
        skipped_header = False
        for i, line in enumerate(infile, start=1):
        
            if comment_character is not None and line.startswith(comment_character):
                continue

            if has_header_row and not skipped_header:
                skipped_header = True
                continue

            if filter_field is not None:
                filter_field_value = line[:-1].split(delim)[filter_field]
                if filter_field_value not in filter_set:
                    continue

            self.load_metadata['record_counter'] += 1
            try:
                # TODO we should pass the whole file iterator to csv reader, just need to handle comments and headers
                # CSV Reader expects a list of rows as input and outputs a list of strings.
                # We process one at a time, so we pass "line" in as a list and take the first result.
                reader = csv.reader([line], delimiter=delim)
                split_row = list(reader)[0]
                self.parse_row(split_row, subject_extractor, object_extractor, predicate_extractor, subject_property_extractor, object_property_extractor, edge_property_extractor, exclude_unconnected_nodes)
            except Exception as e:
                self.load_metadata['errors'].append(e.__str__())
                self.load_metadata['skipped_record_counter'] += 1

    def sql_extract(self, cursor, sql_query, subject_extractor, object_extractor, predicate_extractor, subject_property_extractor, object_property_extractor, edge_property_extractor):
        """Read a csv, perform callbacks to retrieve node and edge info per row.
        Assumes that all of the properties extractable for a node occur on the line with the node identifier"""

        cursor.execute(sql_query)
        rows = cursor.fetchall()
        for row in rows:
            self.load_metadata['record_counter'] += 1
            try:
                self.parse_row(row, subject_extractor, object_extractor, predicate_extractor, subject_property_extractor, object_property_extractor, edge_property_extractor)
            except Exception as e:
                self.load_metadata['errors'].append(e.__str__())
                self.load_metadata['skipped_record_counter'] += 1

    def json_extract(self,
                     json_array,
                     subject_extractor,
                     object_extractor,
                     predicate_extractor,
                     subject_property_extractor,
                     object_property_extractor,
                     edge_property_extractor):
        for item in json_array:
            self.load_metadata['record_counter'] += 1
            try:
                self.parse_row(item, subject_extractor, object_extractor, predicate_extractor, subject_property_extractor, object_property_extractor, edge_property_extractor)
            except Exception as e:
                self.load_metadata['errors'].append(e.__str__())
                self.load_metadata['skipped_record_counter'] += 1
                return

    def process_tuple(self,
                      subject_id: str, object_id: Optional[str], predicate: Optional[str],
                      subjectprops: Optional[Mapping[str, str]], objectprops: Optional[Mapping[str, str]],
                      edgeprops: Optional[Mapping[str, str]], exclude_unconnected_nodes=False):

        if exclude_unconnected_nodes and not predicate:
            return
        self.load_metadata['record_counter'] += 1
        try:
            subjectprops = subjectprops or {}
            objectprops = objectprops or {}
            edgeprops = edgeprops or {}

            # if we  haven't seen the subject before, add it to nodes
            if subject_id and subject_id not in self.node_ids:
                subject_name = subjectprops.pop('name', '')
                subject_categories = subjectprops.pop('categories', None)
                subject_node = kgxnode(subject_id, name=subject_name, categories=subject_categories, nodeprops=subjectprops)
                if self.file_writer:
                    self.file_writer.write_kgx_node(subject_node)
                else:
                    self.nodes.append(subject_node)
                    self.node_ids.add(subject_id)

            if subject_id and object_id and predicate:
                primary_knowledge_source = edgeprops.pop(PRIMARY_KNOWLEDGE_SOURCE, None)
                aggregator_knowledge_sources = edgeprops.pop(AGGREGATOR_KNOWLEDGE_SOURCES, None)
                edge = kgxedge(subject_id,
                            object_id,
                            predicate=predicate,
                            primary_knowledge_source=primary_knowledge_source,
                            aggregator_knowledge_sources=aggregator_knowledge_sources,
                            edgeprops=edgeprops)
                if self.file_writer:
                    self.file_writer.write_kgx_edge(edge)
                else:
                    self.edges.append(edge)
        except Exception as e:
            self.load_metadata['errors'].append(e.__str__())
            self.load_metadata['skipped_record_counter'] += 1

    def parse_row(self,
                  row,
                  subject_extractor,
                  object_extractor,
                  predicate_extractor,
                  subject_property_extractor,
                  object_property_extractor,
                  edge_property_extractor,
                  exclude_unconnected_nodes=False):
        # pull the information out of the edge
        predicate = predicate_extractor(row) if predicate_extractor is not None else None
        if exclude_unconnected_nodes and predicate is None:
            return
        subject_id = subject_extractor(row)
        object_id = object_extractor(row) if object_extractor is not None else None
        subjectprops = subject_property_extractor(row) if subject_property_extractor is not None else {}
        objectprops = object_property_extractor(row) if object_property_extractor is not None else {}
        edgeprops = edge_property_extractor(row) if edge_property_extractor is not None else {}

        self.process_tuple(subject_id, object_id, predicate, subjectprops, objectprops, edgeprops, exclude_unconnected_nodes)

    def get_node_ids(self):
        if self.file_writer:
            return self.file_writer.written_nodes
        else:
            return self.node_ids

