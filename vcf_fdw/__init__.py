#!/usr/bin/env python
#
# -*- coding: utf-8 -*-
#
# vcf_fdw - a multicorn based foreign data wrapper for VCF 4.0 / 4.1 files using PyVCF
#
# Copyright (c) 2016 Ernst-Georg Schmid
#
# Permission to use, copy, modify, and distribute this software and its documentation for any purpose, without fee, and without a written agreement is hereby granted, provided that the above copyright notice and this paragraph and the following two paragraphs appear in all copies.
#
# IN NO EVENT SHALL Ernst-Georg Schmid BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF Ernst-Georg Schmid HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Ernst-Georg Schmid SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
# THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND Ernst-Georg Schmid HAS NO OBLIGATIONS TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.


from multicorn import ForeignDataWrapper, ANY, TableDefinition, ColumnDefinition
from multicorn.utils import log_to_postgres, ERROR, DEBUG
import vcf
import os


__author__ = "Ernst-Georg Schmid"
__copyright__ = "Copyright (c) 2016 Ernst-Georg Schmid"
__license__ = "PostgreSQL"
__version__ = "1.0"


class VCFForeignDataWrapper(ForeignDataWrapper):
    """A multicorn based foreign data wrapper for VCF 4.x files using pycvf"""

    def __init__(self, options, columns):
        super(VCFForeignDataWrapper, self).__init__(options, columns)
        self.basedir = options.get('basedir', '')
        self.suffix = options.get('suffix', '.vcf.gz')
        self.species = options.get('species', '')

    @classmethod
    def import_schema(
            self,
            schema,
            srv_options,
            options,
            restriction_type,
            restricts):
        """Support for IMPORT FOREIGN SCHEMA"""
        ftable = TableDefinition(
            options['species'], schema=None, columns=[
                ColumnDefinition(
                    "chrom", type_name="text"), ColumnDefinition(
                    "pos", type_name="integer"), ColumnDefinition(
                    "id", type_name="text"), ColumnDefinition(
                        "ref", type_name="text"), ColumnDefinition(
                            "alt", type_name="text[]"), ColumnDefinition(
                                "qual", type_name="real"), ColumnDefinition(
                                    "heterozygosity", type_name="real"), ColumnDefinition(
                                        "sample", type_name="text"), ColumnDefinition(
                                            "species", type_name="text"), ColumnDefinition(
                                                "info", type_name="text"), ColumnDefinition(
                                                    "depth", type_name="integer"),
                ColumnDefinition(
                    "genotype", type_name="text"), ColumnDefinition(
                    "filter", type_name="text"), ColumnDefinition(
                    "issnp", type_name="boolean"), ColumnDefinition(
                    "issv", type_name="boolean"), ColumnDefinition(
                    "isindel", type_name="boolean"), ColumnDefinition(
                    "ismonomorphic", type_name="boolean"), ColumnDefinition(
                    "isdeletion", type_name="boolean"), ColumnDefinition(
                    "issvprecise", type_name="boolean"), ColumnDefinition(
                    "istransition", type_name="boolean"), ColumnDefinition(
                    "source", type_name="text")])

        ftable.options = options
        ftable.options['basedir'] = schema

        return [ftable]
        # e.g. import foreign schema "/home/ergo/export" from server
        # multicorn_vcf into vcf options (species 'human', suffix
        # '_vuid.vcf.gz')

    def get_rel_size(self, quals, columns):
        """Inform the planner about the cost of a query"""
        width = len(columns) * 16
        nb_rows = 10000

        return (nb_rows, width)

    def get_path_keys(self):
        """Helps the planner by supplying a list of list of access keys,
        as well as a row estimate for each one."""
        return [(('sample',), 10000000), (('chrom',), 10000000), (('pos',), 100000), (('id',), 10000), ((
            'chrom', 'pos'), 1000), (('chrom', 'pos', 'id'), 100), (('sample', 'chrom', 'pos', 'id'), 1)]

    def execute(self, quals, columns, sortkeys=None):
        """Execute a query"""
        if __name__ == 'vcf_fdw':
            if not quals:
                log_to_postgres(
                    'Query not supported: NO_QUALS',
                    ERROR,
                    'You may add some qualifiers.')

            sample = []
            chrom = []
            positions = []

            if quals:
                for qual in quals:
                    if qual.field_name.lower() == 'sample':
                        if qual.operator == '=':
                            sample.append(qual.value)
                        elif qual.list_any_or_all is ANY:
                            if qual.operator[0] == '=':
                                sample.extend(qual.value)

                    elif qual.field_name.lower() == 'chrom':
                        if qual.operator == '=':
                            chrom.append(qual.value)
                        elif qual.list_any_or_all is ANY:
                            if qual.operator[0] == '=':
                                chrom.extend(qual.value)

                    elif qual.field_name.lower() == 'pos' and qual.operator == '<=':
                        positions = [(0, qual.value)]
                    elif qual.field_name.lower() == 'pos' and qual.operator == '<':
                        positions = [(0, (qual.value - 1))]
                    elif qual.field_name.lower() == 'pos' and qual.operator == '>=':
                        positions = [(qual.value, None)]
                    elif qual.field_name.lower() == 'pos' and qual.operator == '>':
                        positions = [((qual.value + 1), None)]
                    elif qual.field_name.lower() == 'pos':
                        if qual.operator == '=':
                            positions = [((qual.value - 1), qual.value)]
                        elif qual.list_any_or_all is ANY:
                            if qual.operator[0] == '=':
                                for qv in qual.value:
                                    positions.append((qv - 1, qv))

            if positions:
                positions = deduplicate(positions)

            if not sample:
                for e in os.listdir(os.path.join(self.basedir, self.species)):
                    if e.endswith(self.suffix):
                        p = e.find(self.suffix)
                        if p != -1:
                            sample.append(e[:p])
            else:
                sample = deduplicate(sample)

            for sn in sample:
                filename = os.path.join(
                    self.basedir, self.species, sn + self.suffix)
                try:
                    vcf_reader = vcf.Reader(filename=filename, compressed=True)
                except IOError:
                    log_to_postgres(
                        'Data file {0} for sample {1} not found.'.format(
                            filename, sn), DEBUG, 'Check filesystem.')
                    continue

                fetch = vcf_reader.fetch

                if not chrom:
                    chrom = []
                    for ct in vcf_reader.contigs.iterkeys():
                        chrom.append(ct)

                chrom = to_ascii(deduplicate(chrom))

                if not positions:
                    positions = [(None, None)]

                for c in chrom:
                    for p in positions:
                        try:
                            records = fetch(c, p[0], p[1])
                        except ValueError:
                            log_to_postgres(
                                'Contig {0} not in sample {1}.'.format(
                                    c, sn), DEBUG, 'Check VCF.')
                            continue
                        else:
                            for record in records:
                                yield transcode(
                                    columns,
                                    filename,
                                    record,
                                    sn,
                                    self.species)


def transcode(columns, filename, record, sample, species):
    """Make a PostgreSQL record from a pyvcf VCF record"""

    line = {}

    for column_name in columns:
        if column_name.lower() == 'chrom':
            line[column_name] = '{0}'.format(record.CHROM)
        elif column_name.lower() == 'pos':
            line[column_name] = '{0}'.format(record.POS)
        elif column_name.lower() == 'id':
            #_id = '%s' % (record.ID)
            line[column_name] = '{0}'.format(record.ID)
        elif column_name.lower() == 'ref':
            line[column_name] = '{0}'.format(record.REF)
        elif column_name.lower() == 'alt':
            alt = '{0}'.format(record.ALT)
            line[column_name] = alt.replace(
                '[', '{').replace(']', '}').replace("'", '"')
        elif column_name.lower() == 'qual':
            line[column_name] = '{0}'.format(record.QUAL)
        elif column_name.lower() == 'heterozygosity':
            #_het = '%f' % (record.heterozygosity)
            line[column_name] = '{0}'.format(record.heterozygosity)
        elif column_name.lower() == 'filter':
            line[column_name] = '{0}'.format(
                (record.FILTER if record.FILTER else "['PASS']"))
        elif column_name.lower() == 'issnp':
            line[column_name] = '{0}'.format(('t' if record.is_snp else 'f'))
        elif column_name.lower() == 'issv':
            line[column_name] = '{0}'.format(('t' if record.is_sv else 'f'))
        elif column_name.lower() == 'isindel':
            line[column_name] = '{0}'.format(('t' if record.is_indel else 'f'))
        elif column_name.lower() == 'ismonomorphic':
            line[column_name] = '{0}'.format(
                ('t' if record.is_monomorphic else 'f'))
        elif column_name.lower() == 'isdeletion':
            line[column_name] = '{0}'.format(
                ('t' if record.is_deletion else 'f'))
        elif column_name.lower() == 'issvprecise':
            line[column_name] = '{0}'.format(
                ('t' if record.is_sv_precise else 'f'))
        elif column_name.lower() == 'istransition':
            line[column_name] = '{0}'.format(
                ('t' if record.is_transition else 'f'))
        elif column_name.lower() == 'source':
            line[column_name] = '{0}'.format(filename)
        elif column_name.lower() == 'sample':
            line[column_name] = '{0}'.format(sample)
        elif column_name.lower() == 'species':
            line[column_name] = '{0}'.format(species)
        elif column_name.lower() == 'info':
            line[column_name] = '{0}'.format(record.INFO)
        elif column_name.lower() == 'depth':
            line[column_name] = '{0}'.format(record.INFO['DP'])
        elif column_name.lower() == 'genotype':
            try:
                call = record.genotype(sample)
            except KeyError:
                line[column_name] = '{0}'.format('N/A')
            else:
                line[column_name] = '{0}'.format(call['GT'])

    return line


def deduplicate(vals):
    """Deduplicate values in a list"""
    return [v for v in set(vals)]


def to_ascii(vals):
    """Convert values in a list to ASCII"""
    return [v.encode('ascii', 'ignore') for v in vals]
