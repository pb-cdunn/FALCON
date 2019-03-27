version 1.0
task task__0_rawreads__build {
  input {
    File config
    Array[File] fasta_array
    String topdir = "../.."
    String pypeflow_nproc = "1"
    String pypeflow_mb = "4000"
  }
  command <<<
    files_fn="~{write_lines(fasta_array)}"
    python -m falcon_kit.mains.dazzler --config-fn=~{config} --db-fn=./raw_reads.db  build --input-fofn-fn=${files_fn} --length-cutoff-fn=./length_cutoff
    # TODO: Verify that db exists.
    #ln -sf ./length_cutoff length_cutoff
  >>>
  output {
    File dazzler_db = "raw_reads.db"
    File dazzler_idx = ".raw_reads.idx"
    File dazzler_bps = ".raw_reads.bps"
    File dazzler_dust_anno = ".raw_reads.dust.anno"
    File dazzler_dust_data = ".raw_reads.dust.data"
    File length_cutoff = "length_cutoff"
  }
}
task task__0_rawreads__tan_split {
  input {
    File config
    File dazzler_db = "raw_reads.db"
    File dazzler_idx = ".raw_reads.idx"
    File dazzler_bps = ".raw_reads.bps"
    #File dazzler_dust_anno = ".raw_reads.dust.anno"
    #File dazzler_dust_data = ".raw_reads.dust.data"
    #File length_cutoff
    #String pypeflow_nproc = "1"
    #String wildcards = "dal0_id"
  }
  command <<<
    python -m falcon_kit.mains.dazzler --config=~{config} --db=~{dazzler_db}                            tan-split --split=./tan-uows.json --bash-template=./bash_template.sh
  >>>
  output {
    File uows = 'all-units-of-work.tar'
  }
}
task task__0_rawreads__tan_scatter {
  input {
    File uows
  }
  command <<<
    python -m falcon_kit.mains.generic_scatter_uows_tar --all=~{uows} --nchunks=16 --pattern='./some-units-of-work.%.tar'
  >>>
  output {
    Array[File] chunks = glob("some-units-of-work.*.tar")
  }
}
task task__0_rawreads__tan_apply {
  input {
    File chunk
    File config
    File dazzler_db
    File dazzler_idx
    File dazzler_bps
    #File dazzler_dust_anno
    #File dazzler_dust_data
    Int nproc = 4
  }
  command <<<
    set -vex
    python -m falcon_kit.mains.cromwell_symlink ~{config} ~{dazzler_db} ~{dazzler_idx} ~{dazzler_bps}
    python -m falcon_kit.mains.cromwell_run_uows_tar --nproc=~{nproc} --uows-tar-fn=~{chunk} --tool=datander
    # glob seems to fail on dotfiles, so we rename them first.
    ls -la uow-*/.raw_reads.*.tan.anno
    ls -la uow-*/.raw_reads.*.tan.data
    python -m falcon_kit.mains.cromwell_undot --pattern 'uow-*/.*.*.tan.anno' --prefix dot
    python -m falcon_kit.mains.cromwell_undot --pattern 'uow-*/.*.*.tan.data' --prefix dot
    #find . -name '.*.*.tan.anno' | xargs -I XXX mv -sf XXX dotXXX
    #find . -name '.*.*.tan.data' | xargs -I XXX ln -sf XXX dotXXX
    #ls -larth .
    #Catrack -vdf raw_reads tan
  >>>
  output {
    Array[File] dazzler_tan_anno_array = glob("uow-*/dot.*.*.tan.anno")
    Array[File] dazzler_tan_data_array = glob("uow-*/dot.*.*.tan.data")
  }
}
task task__0_rawreads__tan_combine {
  input {
    File config
    File dazzler_db
    File dazzler_idx
    File dazzler_bps
    Array[Array[File]] tan_anno_array_of_array
    Array[Array[File]] tan_data_array_of_array
  }
  Array[File] tan_anno_array = flatten(tan_anno_array_of_array)
  Array[File] tan_data_array = flatten(tan_data_array_of_array)
  # There are problems with write_json: https://github.com/broadinstitute/cromwell/issues/4625

  command <<<
    set -vex
    anno_files_fn="~{write_lines(tan_anno_array)}"
    data_files_fn="~{write_lines(tan_data_array)}"
    cat ${anno_files_fn}

    python -m falcon_kit.mains.dazzler --config=~{config} --db=~{dazzler_db}  track-combine --track tan --anno ${anno_files_fn} --data ${data_files_fn}
    # Note: track-combine, not tan-combine, which was for the UOW workflow
  >>>
  output {
    File dazzler_tan_anno = ".raw_reads.tan.anno"
    File dazzler_tan_data = ".raw_reads.tan.data"
  }
}
task task__0_rawreads__daligner_split {
  input {
    File config
    File dazzler_db = "raw_reads.db"
    File dazzler_idx = ".raw_reads.idx"
    File dazzler_bps = ".raw_reads.bps"
    File dazzler_dust_anno = ".raw_reads.dust.anno"
    File dazzler_dust_data = ".raw_reads.dust.data"
    #File dazzler_tan_anno = ".raw_reads.tan.anno"
    #File dazzler_tan_data = ".raw_reads.tan.data"
    File length_cutoff
    String pypeflow_nproc = "1"
    String wildcards = "dal0_id"
  }
  command <<<
    python -m falcon_kit.mains.dazzler --config=~{config} --db=~{dazzler_db} --nproc=~{pypeflow_nproc}  daligner-split --wildcards=~{wildcards} --length-cutoff-fn=~{length_cutoff} --split-fn=./all-units-of-work.json --bash-template-fn=./daligner_bash_template.sh
    #python -m falcon_kit.mains.generic_tar_uows --all=./all-units-of-work.json --nchunks=16 --pattern=./units-of-work.%.tar
    #python -m falcon_kit.mains.generic_scatter_one_uow --all-uow-list-fn=./all-units-of-work.json --one-uow-list-fn=some-units-of-work.json --split-idx=0
  >>>
  output {
    File uows = 'all-units-of-work.tar'
  }
}
task task__0_rawreads__daligner_scatter {
  input {
    File uows
  }
  command <<<
    python -m falcon_kit.mains.generic_scatter_uows_tar --all=~{uows} --nchunks=16 --pattern='./some-units-of-work.%.tar'
  >>>
  output {
    Array[File] chunks = glob("some-units-of-work.*.tar")
  }
}
task task__0_rawreads__daligner_apply {
  input {
    File chunk
    File config
    File dazzler_db
    File dazzler_idx
    File dazzler_bps
    File dazzler_dust_anno
    File dazzler_dust_data
    File dazzler_tan_anno
    File dazzler_tan_data
    Int nproc = 4
  }
  command <<<
    set -e
    python -m falcon_kit.mains.cromwell_symlink ~{config} ~{dazzler_db} ~{dazzler_idx} ~{dazzler_bps} ~{dazzler_dust_anno} ~{dazzler_dust_data} ~{dazzler_tan_anno} ~{dazzler_tan_data}
    #echo 'python -m falcon_kit.mains.dazzler --config=~{config} --db=~{dazzler_db}  daligner-apply --script=./run_daligner.sh --job-done={output.job_done}' > ./bash_template.sh

    python -m falcon_kit.mains.cromwell_run_uows_tar --nproc=~{nproc} --uows-tar-fn=~{chunk} --tool=daligner
  >>>
  output {
    Array[File] result = glob("uow-*/*.las")
  }
}
task task__0_rawreads__daligner_las_merge {
  input {
    Array[Array[File]] las_array_of_array
    File config
  }
  Array[File] las_fns = flatten(las_array_of_array)
  # There are problems with write_json: https://github.com/broadinstitute/cromwell/issues/4625

  command <<<
    set -vex
    files_fn="~{write_lines(las_fns)}"
    echo ${files_fn}
    python -m falcon_kit.mains.cromwell_write_json --lines-fn=${files_fn} --json-fn=./gathered-las.json

    python -m falcon_kit.mains.cromwell_symlink ~{config}

    python -m falcon_kit.mains.dazzler --config=~{basename(config)}                  merge-split --db-prefix=raw_reads --las-paths=./gathered-las.json --wildcards=mer0_id --split-fn=all-units-of-work.json --bash-template-fn=las-merge-bash-template.sh

    python -m falcon_kit.mains.generic_run_units_of_work --nproc=1 --bash-template-fn=./las-merge-bash-template.sh --units-of-work-fn=./all-units-of-work.json --results-fn=./merged-las-paths.json
  >>>
  output {
    Array[File] las_array = glob("uow-*/*.las")
  }
}
task task__0_rawreads__cns_apply {
  input {
    #File chunk
    File las
    File config
    File dazzler_db
    File dazzler_idx
    File dazzler_bps
    File length_cutoff
    Int nproc = 8
  }
  command <<<
  set -vex
  python -m falcon_kit.mains.cromwell_symlink ~{config} ~{dazzler_db} ~{dazzler_idx} ~{dazzler_bps} ~{length_cutoff} ~{las}

  python -m falcon_kit.mains.consensus_task --nproc=~{nproc} --las-fn=~{basename(las)} --db-fn=~{basename(dazzler_db)} --length-cutoff-fn=~{basename(length_cutoff)} --config-fn=~{basename(config)} --fasta-fn=consensus.chunk.fasta
  >>>
  output {
    File fasta = "consensus.chunk.fasta"
  }
}
task task__0_rawreads__report {
  input {
    Array[File] fasta_array
    File config
    File dazzler_db
    File dazzler_idx
    File dazzler_bps
    File length_cutoff
  }
  command <<<
  set -vex
  python -m falcon_kit.mains.cromwell_symlink ~{config} ~{dazzler_db} ~{dazzler_idx} ~{dazzler_bps} ~{length_cutoff}

  files_fn="~{write_lines(fasta_array)}"
  ln -sf ${files_fn} input_preads.fofn

  python -m falcon_kit.mains.task_report_pre_assembly --config-fn=~{basename(config)} --length-cutoff-fn=~{basename(length_cutoff)} --raw-reads-db-fn=~{basename(dazzler_db)} --preads-fofn-fn=./input_preads.fofn --pre-assembly-report-fn=./pre_assembly_stats.json
  >>>
  output {
    File stats = "pre_assembly_stats.json"
  }
}
task task__1_preads_ovl__build {
  input {
    File config
    Array[File] fasta_array
    String topdir = "../.."
    String pypeflow_nproc = "1"
    String pypeflow_mb = "4000"
  }
  command <<<
    files_fn="~{write_lines(fasta_array)}"
    python -m falcon_kit.mains.dazzler --config-fn=~{config} --db-fn=./preads.db  build --input-fofn-fn=${files_fn} --length-cutoff-fn=./length_cutoff
    # TODO: Verify that db exists.
    #ln -sf ./length_cutoff length_cutoff
  >>>
  output {
    File dazzler_db = "preads.db"
    File dazzler_idx = ".preads.idx"
    File dazzler_bps = ".preads.bps"
    File dazzler_dust_anno = ".preads.dust.anno"
    File dazzler_dust_data = ".preads.dust.data"
    File length_cutoff = "length_cutoff"
  }
}
task task__1_preads_ovl__daligner_split {
  input {
    File config
    File dazzler_db
    File dazzler_idx
    File dazzler_bps
    File dazzler_dust_anno
    File dazzler_dust_data
    File length_cutoff
    String pypeflow_nproc = "1"
    String wildcards = "dal0_id"
  }
  command <<<
    python -m falcon_kit.mains.dazzler --config=~{config} --db=~{dazzler_db} --nproc=~{pypeflow_nproc}  daligner-split --wildcards=~{wildcards} --length-cutoff-fn=~{length_cutoff} --split-fn=./all-units-of-work.json --bash-template-fn=./daligner_bash_template.sh
    #python -m falcon_kit.mains.generic_tar_uows --all=./all-units-of-work.json --nchunks=16 --pattern=./units-of-work.%.tar
    #python -m falcon_kit.mains.generic_scatter_one_uow --all-uow-list-fn=./all-units-of-work.json --one-uow-list-fn=some-units-of-work.json --split-idx=0
  >>>
  output {
    File uows = 'all-units-of-work.tar'
  }
}
task task__1_preads_ovl__daligner_scatter {
  input {
    File uows
  }
  command <<<
    python -m falcon_kit.mains.generic_scatter_uows_tar --all=~{uows} --nchunks=16 --pattern='./some-units-of-work.%.tar'
  >>>
  output {
    Array[File] chunks = glob("some-units-of-work.*.tar")
  }
}
task task__1_preads_ovl__daligner_apply {
  input {
    File chunk
    File config
    File dazzler_db
    File dazzler_idx
    File dazzler_bps
    File dazzler_dust_anno
    File dazzler_dust_data
    #File dazzler_tan_anno
    #File dazzler_tan_data
    Int nproc = 4
  }
  command <<<
    set -e
    python -m falcon_kit.mains.cromwell_symlink ~{config} ~{dazzler_db} ~{dazzler_idx} ~{dazzler_bps} ~{dazzler_dust_anno} ~{dazzler_dust_data}
    #echo 'python -m falcon_kit.mains.dazzler --config=~{config} --db=~{dazzler_db}  daligner-apply --script=./run_daligner.sh --job-done={output.job_done}' > ./bash_template.sh

    python -m falcon_kit.mains.cromwell_run_uows_tar --nproc=~{nproc} --uows-tar-fn=~{chunk} --tool=daligner
  >>>
  output {
    Array[File] result = glob("uow-*/*.las")
  }
}
task task__1_preads_ovl__daligner_las_merge {
  input {
    Array[Array[File]] las_array_of_array
    File config
  }
  Array[File] las_fns = flatten(las_array_of_array)
  # There are problems with write_json: https://github.com/broadinstitute/cromwell/issues/4625

  command <<<
    set -vex
    files_fn="~{write_lines(las_fns)}"
    echo ${files_fn}
    python -m falcon_kit.mains.cromwell_write_json --lines-fn=${files_fn} --json-fn=./gathered-las.json

    python -m falcon_kit.mains.cromwell_symlink ~{config}

    python -m falcon_kit.mains.dazzler --config=~{basename(config)}                  merge-split --db-prefix=preads --las-paths=./gathered-las.json --wildcards=mer0_id --split-fn=all-units-of-work.json --bash-template-fn=las-merge-bash-template.sh

    python -m falcon_kit.mains.generic_run_units_of_work --nproc=1 --bash-template-fn=./las-merge-bash-template.sh --units-of-work-fn=./all-units-of-work.json --results-fn=./merged-las-paths.json
  >>>
  output {
    Array[File] las_array = glob("uow-*/*.las")
  }
}
task task__1_preads_ovl__db2falcon {
  input {
    File dazzler_db
    File dazzler_idx
    File dazzler_bps
  }
  command <<<
    set -vex
    python -m falcon_kit.mains.cromwell_symlink ~{dazzler_db} ~{dazzler_idx} ~{dazzler_bps}
    time DB2Falcon -U ~{basename(dazzler_db)}
    [ -f ./preads4falcon.fasta ] || exit 1
  >>>
  output {
    File preads4falcon_fasta = './preads4falcon.fasta'
  }
}
task task__2_asm_falcon {
  input {
    File dazzler_db
    File dazzler_idx
    File dazzler_bps
    File preads4falcon_fasta
    Array[File] las_array
    String overlap_filtering_setting
    String length_cutoff_pr
    String fc_ovlp_to_graph_option
  }
  command <<<
    set -vex
    python -m falcon_kit.mains.cromwell_symlink ~{dazzler_db} ~{dazzler_idx} ~{dazzler_bps} ~{preads4falcon_fasta}

    files_fn="~{write_lines(las_array)}"
    python -m falcon_kit.mains.cromwell_write_json --lines-fn=${files_fn} --json-fn=./las_fofn.json

    # Given, las_fofn.json,
    # write preads.ovl:

    #overlap_filtering_setting='--max-diff 10000 --max-cov 100000 --min-cov 1 --min-len 1 --bestn 1000 --n-core 0'
    #length_cutoff_pr='1'

    time python -m falcon_kit.mains.ovlp_filter --db ~{basename(dazzler_db)} --las-fofn ./las_fofn.json ~{overlap_filtering_setting} --min-len ~{length_cutoff_pr} --out-fn preads.ovl

    # Given preads.ovl,
    # write sg_edges_list, c_path, utg_data, ctg_paths.
    #fc_ovlp_to_graph_option='--min-len 1'
    time python -m falcon_kit.mains.ovlp_to_graph ~{fc_ovlp_to_graph_option} --overlap-file preads.ovl >| fc_ovlp_to_graph.log

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
  >>>
  output {
    File p_ctg_fasta = "p_ctg.fa"
    File a_ctg_fasta = "a_ctg.fa"
  }
}
workflow falcon {
  input {
    Array[File] array_of_fasta
    #File ifile_config = "/scratch/cdunn/pbcromwell/General_config.json"
    File ifile_config = "General_config.json"
    Int nproc = 1
  }
  Map[String, String] cfg = read_json(ifile_config)

  call task__0_rawreads__build {
    input:
      topdir         = "../..",
      config         = ifile_config,
      fasta_array    = array_of_fasta,
      pypeflow_nproc = "~{nproc}",
      pypeflow_mb    = "4000",
  }
  call task__0_rawreads__tan_split {
    input:
      config            = ifile_config,
      dazzler_db        = task__0_rawreads__build.dazzler_db,
      dazzler_idx       = task__0_rawreads__build.dazzler_idx,
      dazzler_bps       = task__0_rawreads__build.dazzler_bps,
      #dazzler_dust_anno = task__0_rawreads__build.dazzler_dust_anno,
      #dazzler_dust_data = task__0_rawreads__build.dazzler_dust_data,
  }
  call task__0_rawreads__tan_scatter {
    input: uows = task__0_rawreads__tan_split.uows
  }
  scatter (chunk in task__0_rawreads__tan_scatter.chunks) {
    call task__0_rawreads__tan_apply {
      input:
        chunk             = chunk,
        config            = ifile_config,
        dazzler_db        = task__0_rawreads__build.dazzler_db,
        dazzler_idx       = task__0_rawreads__build.dazzler_idx,
        dazzler_bps       = task__0_rawreads__build.dazzler_bps,
        #dazzler_dust_anno = task__0_rawreads__build.dazzler_dust_anno,
        #dazzler_dust_data = task__0_rawreads__build.dazzler_dust_data,
    }
  }
  call task__0_rawreads__tan_combine {
    input:
      config            = ifile_config,
      dazzler_db        = task__0_rawreads__build.dazzler_db,
      dazzler_idx       = task__0_rawreads__build.dazzler_idx,
      dazzler_bps       = task__0_rawreads__build.dazzler_bps,
      tan_anno_array_of_array = task__0_rawreads__tan_apply.dazzler_tan_anno_array,
      tan_data_array_of_array = task__0_rawreads__tan_apply.dazzler_tan_data_array,
  }
  call task__0_rawreads__daligner_split {
    input:
      config            = ifile_config,
      dazzler_db        = task__0_rawreads__build.dazzler_db,
      dazzler_idx       = task__0_rawreads__build.dazzler_idx,
      dazzler_bps       = task__0_rawreads__build.dazzler_bps,
      dazzler_dust_anno = task__0_rawreads__build.dazzler_dust_anno,
      dazzler_dust_data = task__0_rawreads__build.dazzler_dust_data,
      length_cutoff     = task__0_rawreads__build.length_cutoff
  }
  call task__0_rawreads__daligner_scatter {
    input: uows = task__0_rawreads__daligner_split.uows
  }
  scatter (chunk in task__0_rawreads__daligner_scatter.chunks) {
    call task__0_rawreads__daligner_apply {
      input:
        chunk             = chunk,
        config            = ifile_config,
        dazzler_db        = task__0_rawreads__build.dazzler_db,
        dazzler_idx       = task__0_rawreads__build.dazzler_idx,
        dazzler_bps       = task__0_rawreads__build.dazzler_bps,
        dazzler_dust_anno = task__0_rawreads__build.dazzler_dust_anno,
        dazzler_dust_data = task__0_rawreads__build.dazzler_dust_data,
        dazzler_tan_anno = task__0_rawreads__tan_combine.dazzler_tan_anno,
        dazzler_tan_data = task__0_rawreads__tan_combine.dazzler_tan_data,
    }
  }
  call task__0_rawreads__daligner_las_merge {
    input:
      las_array_of_array = task__0_rawreads__daligner_apply.result,
      config             = ifile_config,
  }
  scatter (chunk in task__0_rawreads__daligner_las_merge.las_array) {
    call task__0_rawreads__cns_apply {
      input:
        las           = chunk,
        config        = ifile_config,
        dazzler_db    = task__0_rawreads__build.dazzler_db,
        dazzler_idx   = task__0_rawreads__build.dazzler_idx,
        dazzler_bps   = task__0_rawreads__build.dazzler_bps,
        length_cutoff = task__0_rawreads__build.length_cutoff,
    }
  }
  call task__0_rawreads__report {
    input:
      config        = ifile_config,
      fasta_array   = task__0_rawreads__cns_apply.fasta,
      dazzler_db    = task__0_rawreads__build.dazzler_db,
      dazzler_idx   = task__0_rawreads__build.dazzler_idx,
      dazzler_bps   = task__0_rawreads__build.dazzler_bps,
      length_cutoff = task__0_rawreads__build.length_cutoff,
  }
  call task__1_preads_ovl__build {
    input:
      topdir        = "../..",
      config        = ifile_config,
      fasta_array   = task__0_rawreads__cns_apply.fasta,
      pypeflow_nproc= "~{nproc}",
      pypeflow_mb   = "4000",
  }
  call task__1_preads_ovl__daligner_split {
    input:
      config=ifile_config,
      dazzler_db=task__1_preads_ovl__build.dazzler_db,
      dazzler_idx=task__1_preads_ovl__build.dazzler_idx,
      dazzler_bps=task__1_preads_ovl__build.dazzler_bps,
      dazzler_dust_anno=task__1_preads_ovl__build.dazzler_dust_anno,
      dazzler_dust_data=task__1_preads_ovl__build.dazzler_dust_data,
      length_cutoff=task__1_preads_ovl__build.length_cutoff
  }
  call task__1_preads_ovl__daligner_scatter {
    input: uows = task__1_preads_ovl__daligner_split.uows
  }
  scatter (chunk in task__1_preads_ovl__daligner_scatter.chunks) {
    call task__1_preads_ovl__daligner_apply {
      input:
        chunk = chunk,
        config = ifile_config,
        dazzler_db = task__1_preads_ovl__build.dazzler_db,
        dazzler_idx= task__1_preads_ovl__build.dazzler_idx,
        dazzler_bps= task__1_preads_ovl__build.dazzler_bps,
        dazzler_dust_anno=task__1_preads_ovl__build.dazzler_dust_anno,
        dazzler_dust_data=task__1_preads_ovl__build.dazzler_dust_data,
    }
  }
  call task__1_preads_ovl__daligner_las_merge {
    input:
      las_array_of_array = task__1_preads_ovl__daligner_apply.result,
      config             = ifile_config,
  }
  call task__1_preads_ovl__db2falcon {
    input:
      dazzler_db = task__1_preads_ovl__build.dazzler_db,
      dazzler_idx= task__1_preads_ovl__build.dazzler_idx,
      dazzler_bps= task__1_preads_ovl__build.dazzler_bps,
  }
  call task__2_asm_falcon {
    input:
      dazzler_db          = task__1_preads_ovl__build.dazzler_db,
      dazzler_idx         = task__1_preads_ovl__build.dazzler_idx,
      dazzler_bps         = task__1_preads_ovl__build.dazzler_bps,
      las_array           = task__1_preads_ovl__daligner_las_merge.las_array,
      preads4falcon_fasta = task__1_preads_ovl__db2falcon.preads4falcon_fasta,
      overlap_filtering_setting = cfg["overlap_filtering_setting"],
      length_cutoff_pr          = cfg["length_cutoff_pr"],
      fc_ovlp_to_graph_option   = cfg["fc_ovlp_to_graph_option"],
  }
  output {
    File ofile_pre_assembly_report = task__0_rawreads__report.stats
    File ofile_preads4falcon_fasta = task__1_preads_ovl__db2falcon.preads4falcon_fasta
    File ofile_p_ctg_fasta         = task__2_asm_falcon.p_ctg_fasta
    File ofile_a_ctg_fasta         = task__2_asm_falcon.a_ctg_fasta
  }
}
