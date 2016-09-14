import os.path
import re

def get_rule(name):
    return getattr(rules, name)

def is_target_rule(name):
    return not get_rule(name).has_wildcards()

rule all_rulegraphs:
    input:
        expand('rulegraphs/{graphtype}-{target}.{filetype}',
               graphtype=['rulegraph',],
               target=['all', 'all_rnaseq_counts', 'all_salmon', 'all_kallisto', 'all_chipseq_bai', 'all_macs_callpeak', 'all_epic_callpeak',],
               filetype=['pdf', 'png'])

rule all_dags:
    input:
        expand('rulegraphs/{graphtype}-{target}.{filetype}',
               graphtype=['dag',],
               target=['all', 'all_rnaseq_counts', 'all_salmon', 'all_kallisto', 'all_chipseq_bai', 'all_macs_callpeak', 'all_epic_callpeak',],
               filetype=['pdf', 'png'])

rule svg_to_png:
    input: '{filename}.svg'
    output: '{filename}.png'
    shell: '''inkscape {input:q} --export-png={output:q} --export-dpi=300'''

rule svg_to_pdf:
    input: '{filename}.svg'
    output: '{filename}.pdf'
    shell: '''inkscape {input:q} --export-pdf={output:q} --export-dpi=300'''

rule dag_svg:
    input:  'Snakefile', 'rulegraph.Snakefile'
    output: 'rulegraphs/dag-{target}.svg'
    params: target_path=lambda wildcards: re.sub(":", os.path.sep, wildcards.target)
    run:
        if is_target_rule(params.target_path):
            rule = get_rule(params.target_path)
            if len(rule.output):
                real_targets = rule.output
            else:
                real_targets = rule.input
        else:
            real_targets = [params.target_path]
        shell('''
        set -o pipefail;
        snakemake --nolock -f --dag {real_targets:q} | \
        dot -Grankdir=LR -Tsvg > {output:q}
        ''')

rule rulegraph_svg:
    input: 'Snakefile', 'rulegraph.Snakefile'
    output: 'rulegraphs/rulegraph-{target}.svg'
    params: target_path=lambda wildcards: re.sub(":", os.path.sep, wildcards.target)
    run:
        if is_target_rule(params.target_path):
            rule = get_rule(params.target_path)
            if len(rule.output):
                real_targets = rule.output
            else:
                real_targets = rule.input
        else:
            real_targets = [params.target_path]
        shell('''
        set -o pipefail;
        snakemake --nolock -f --rulegraph {real_targets:q} | \
        dot -Grankdir=TB -Tsvg > {output:q}
        ''')
