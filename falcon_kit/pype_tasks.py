from __future__ import absolute_import
from __future__ import print_function


from future.utils import viewitems
from future.utils import itervalues
# PypeTask functions now need to be module-level.
from . import run_support as support
from . import bash  # for scattering
# from pypeflow.simple_pwatcher_bridge import fn # not really needed
import collections
import json
import logging
import os.path
LOG = logging.getLogger(__name__)


TASK_BAM2DEXTA_SPLIT_SCRIPT = """\
python -m falcon_kit.mains.bam2dexta --config-fn={input.config} split  --wildcards={params.wildcards} --bam={input.bam} --split-fn={output.split} --bash-template-fn={output.bash_template}
"""
TASK_BAM2DEXTA_APPLY_SCRIPT = """\
python -m falcon_kit.mains.bam2dexta --config-fn={input.config} apply  --bam-fn={input.bam} --dexta-fn={output.dexta}
"""
TASK_BAM2DEXTA_COMBINE_SCRIPT = """\
python -m falcon_kit.mains.bam2dexta --config-fn={input.config} combine  --gathered-fn={input.gathered} --dexta-fofn-fn={output.fofn}
"""
TASK_CONSENSUS_SPLIT_SCRIPT = """\
python -m falcon_kit.mains.consensus_split --wildcards={params.wildcards} --p-id2las-fn={input.p_id2las} --db-fn={input.raw_reads_db} --length-cutoff-fn={input.length_cutoff} --config-fn={input.config} --split-fn={output.split} --bash-template-fn={output.bash_template}
"""
TASK_CONSENSUS_TASK_SCRIPT = """\
python -m falcon_kit.mains.consensus_task --nproc={params.pypeflow_nproc} --las-fn={input.las} --db-fn={input.db} --length-cutoff-fn={input.length_cutoff} --config-fn={input.config} --fasta-fn={output.fasta}
"""
TASK_CONSENSUS_GATHER_SCRIPT = """\
python -m falcon_kit.mains.consensus_gather_fasta_fofn --gathered-fn={input.gathered} --preads-fofn-fn={output.preads_fofn}
"""
TASK_REPORT_PRE_ASSEMBLY_SCRIPT = """\
python -m falcon_kit.mains.task_report_pre_assembly --config-fn={input.config} --length-cutoff-fn={input.length_cutoff} --raw-reads-db-fn={input.raw_reads_db} --preads-fofn-fn={input.preads_fofn} --pre-assembly-report-fn={output.pre_assembly_report}
"""

TASK_DB_BUILD_SCRIPT = """\
python -m falcon_kit.mains.dazzler --config-fn={input.config} --db-fn={output.db}  build --input-fofn-fn={input.input_fofn} --length-cutoff-fn={output.length_cutoff}
# TODO: Verify that db exists.
#ln -sf {output.length_cutoff} length_cutoff
"""
TASK_DB_TAN_SPLIT_SCRIPT = """\
python -m falcon_kit.mains.dazzler --config={input.config} --db={input.db}  tan-split --split={output.split} --bash-template={output.bash_template}
"""
TASK_DB_TAN_APPLY_SCRIPT = """\
python -m falcon_kit.mains.dazzler --config={input.config} --db={input.db}  tan-apply --script={input.script} --job-done={output.job_done}
"""
TASK_DB_TAN_COMBINE_SCRIPT = """\
python -m falcon_kit.mains.dazzler --config={input.config} --db={input.db}  tan-combine --gathered={input.gathered} --new-db={output.new_db}
"""
TASK_DB_REP_SPLIT_SCRIPT = """\
python -m falcon_kit.mains.dazzler --config={input.config} --db={input.db}  rep-split --las-paths-fn={input.las_paths} --wildcards={params.wildcards} -g{params.group_size} -c{params.coverage_limit} --split={output.split} --bash-template={output.bash_template}
"""
TASK_DB_REP_APPLY_SCRIPT = """\
python -m falcon_kit.mains.dazzler --config={input.config} --db={input.db}  rep-apply --script={input.script} --job-done={output.job_done}
"""
TASK_DB_REP_COMBINE_SCRIPT = """\
python -m falcon_kit.mains.dazzler --config={input.config} --db={input.db}  rep-combine -g{params.group_size} --gathered={input.gathered} --new-db={output.new_db}
"""
TASK_DB_REP_DALIGNER_SPLIT_SCRIPT = """\
python -m falcon_kit.mains.dazzler --config={input.config} --db={input.db} --nproc={params.pypeflow_nproc}  rep-daligner-split --wildcards={params.wildcards} --group-size={params.group_size} --coverage-limit={params.coverage_limit} --split-fn={output.split} --bash-template-fn={output.bash_template}
"""
TASK_DB_DALIGNER_SPLIT_SCRIPT = """\
python -m falcon_kit.mains.dazzler --config={input.config} --db={input.db} --nproc={params.pypeflow_nproc}  daligner-split --wildcards={params.wildcards} --length-cutoff-fn={input.length_cutoff} --split-fn={output.split} --bash-template-fn={output.bash_template}
"""
TASK_DB_DALIGNER_APPLY_SCRIPT = """\
python -m falcon_kit.mains.dazzler --config={input.config} --db={input.db}  daligner-apply --script={input.script} --job-done={output.job_done}
"""
TASK_DB_DALIGNER_COMBINE_SCRIPT = """\
python -m falcon_kit.mains.dazzler --config={input.config} --db={input.db}  daligner-combine --gathered={input.gathered} --las-paths-fn={output.las_paths}
"""
TASK_DB_LAMERGE_SPLIT_SCRIPT = """\
python -m falcon_kit.mains.dazzler --config={input.config}                  merge-split --db-prefix={params.db_prefix} --las-paths={input.las_paths} --wildcards={params.wildcards} --split-fn={output.split} --bash-template-fn={output.bash_template}
"""
TASK_DB_LAMERGE_APPLY_SCRIPT = """\
python -m falcon_kit.mains.dazzler --config={input.config}                  merge-apply --las-paths={input.las_paths} --las-fn={output.las_fn}
"""
TASK_DB_LAMERGE_COMBINE_SCRIPT = """\
python -m falcon_kit.mains.dazzler --config={input.config}                  merge-combine --gathered={input.gathered} --las-paths-fn={output.las_paths} --block2las-fn={output.block2las}
"""

TASK_DUMP_RAWREAD_IDS_SCRIPT = """\
DBshow -n {input.rawread_db} | tr -d '>' | LD_LIBRARY_PATH= awk '{{print $1}}' > {output.rawread_id_file}
"""
TASK_DUMP_PREAD_IDS_SCRIPT = """\
DBshow -n {input.pread_db} | tr -d '>' | LD_LIBRARY_PATH= awk '{{print $1}}' > {output.pread_id_file}
"""
TASK_GENERATE_READ_TO_CTG_MAP_SCRIPT = """\
python -m falcon_kit.mains.generate_read_to_ctg_map --rawread-id={input.rawread_id_file} --pread-id={input.pread_id_file} --sg-edges-list={input.sg_edges_list} --utg-data={input.utg_data} --ctg-paths={input.ctg_paths} --output={output.read_to_contig_map}
"""
TASK_RUN_DB_TO_FALCON_SCRIPT = """\
# Given preads.db,
# write preads4falcon.fasta (implicitly) in CWD.
time DB2Falcon -U {input.preads_db}
[ -f {output.preads4falcon} ] || exit 1
touch {output.job_done}
"""
TASK_RUN_FALCON_ASM_SCRIPT = """\
# Given, las_fofn.json,
# write preads.ovl:

# mobs uses binwrappers, so it does not see our "entry-points".
# So, after dropping "src/py_scripts/*.py", we can call these via python -m:

time python -m falcon_kit.mains.ovlp_filter --db {input.db_file} --las-fofn {input.las_fofn} {params.overlap_filtering_setting} --min-len {params.length_cutoff_pr} --out-fn preads.ovl

ln -sf {input.preads4falcon_fasta} ./preads4falcon.fasta

# Given preads.ovl,
# write sg_edges_list, c_path, utg_data, ctg_paths.
time python -m falcon_kit.mains.ovlp_to_graph {params.fc_ovlp_to_graph_option} --overlap-file preads.ovl >| fc_ovlp_to_graph.log

# Given sg_edges_list, utg_data, ctg_paths, preads4falcon.fasta,
# write p_ctg.fa and a_ctg_all.fa,
# plus a_ctg_base.fa, p_ctg_tiling_path, a_ctg_tiling_path, a_ctg_base_tiling_path:
time python -m falcon_kit.mains.graph_to_contig

# Given a_ctg_all.fa, write a_ctg.fa:
time python -m falcon_kit.mains.dedup_a_tigs >| a_ctg.fa

# Given a_ctg.fa and a_ctg_all_tiling_path, write a_ctg_tiling_path:
time python -m falcon_kit.mains.dedup_a_tp >| a_ctg_tiling_path

# Collect all info needed to format the GFA-1 and GFA-2 representations of
# the assembly graphs.
time python -m falcon_kit.mains.collect_pread_gfa >| asm.gfa.json
time python -m falcon_kit.mains.collect_pread_gfa --add-string-graph >| sg.gfa.json
time python -m falcon_kit.mains.collect_contig_gfa >| contig.gfa.json

# Output the assembly pread graph.
time python -m falcon_kit.mains.gen_gfa_v1 asm.gfa.json >| asm.gfa
time python -m falcon_kit.mains.gen_gfa_v2 asm.gfa.json >| asm.gfa2

# Output the string graph.
time python -m falcon_kit.mains.gen_gfa_v1 sg.gfa.json >| sg.gfa
time python -m falcon_kit.mains.gen_gfa_v2 sg.gfa.json >| sg.gfa2

# Output the contig graph with associate contigs attached to each primary contig.
time python -m falcon_kit.mains.gen_gfa_v2 contig.gfa.json >| contig.gfa2

#rm -f ./preads4falcon.fasta

touch {output.falcon_asm_done}
"""


def fn(p): return p


def system(call, check=False):
    LOG.debug('$(%s)' % repr(call))
    rc = os.system(call)
    msg = 'Call %r returned %d.' % (call, rc)
    if rc:
        LOG.warning(msg)
        if check:
            raise Exception(msg)
    else:
        LOG.debug(msg)
    return rc


def task_dump_rawread_ids(self):
    rawread_db = fn(self.rawread_db)
    rawread_id_file = fn(self.rawread_id_file)
    input = object()
    input.rawread_db = rawread_db
    output = object()
    output.rawread_id_file = rawread_id_file
    system(TASK_DUMP_RAWREAD_IDS_SCRIPT.format(**locals()))


def task_dump_pread_ids(self):
    pread_db = fn(self.pread_db)
    pread_id_file = fn(self.pread_id_file)
    input = object()
    input.pread_db = pread_db
    output = object()
    output.pread_id_file = pread_id_file
    system(TASK_DUMP_PREAD_IDS_SCRIPT.format(**locals()))


def task_generate_read_to_ctg_map(self):
    input = object()
    input.rawread_id_file = fn(self.rawread_id_file)
    input.pread_id_file = fn(self.pread_id_file)
    input.sg_edges_list = fn(self.sg_edges_list)
    input.utg_data = fn(self.utg_data)
    input.ctg_paths = fn(self.ctg_paths)
    output = object()
    output.read_to_contig_map = fn(self.read_to_contig_map)
    system(TASK_GENERATE_READ_TO_CTG_MAP_SCRIPT.format(**locals()))
