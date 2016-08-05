rule rulegraphs:
    input: 'rulegraph-all.pdf', 'dag-all.pdf'

rule svg_to_png:
    input: '{filename}.svg'
    output: '{filename}.png'
    # shell: '''
    # rsvg-convert -u -d 180 -p 180 -f png {input:q} -o {output:q} &&
    #   [ -s {output:q} ] || {{
    #     echo >&2 "ERROR: Output has zero size";
    #     rm -f {output:q};
    #     exit 1;
    #   }}
    # '''
    shell: '''inkscape {input:q} --export-png={output:q} --export-dpi=300'''

rule svg_to_pdf:
    input: '{filename}.svg'
    output: '{filename}.pdf'
    # shell: '''
    # rsvg-convert -u -d 180 -p 180 -f png {input:q} -o {output:q} &&
    #   [ -s {output:q} ] || {{
    #     echo >&2 "ERROR: Output has zero size";
    #     rm -f {output:q};
    #     exit 1;
    #   }}
    # '''
    shell: '''inkscape {input:q} --export-pdf={output:q} --export-dpi=300'''

rule dag_svg:
    input:  'Snakefile', 'rulegraph.Snakefile'
    output: 'dag-{target}.svg'
    shell: '''
    set -o pipefail;
    snakemake --nolock -f --dag {wildcards.target:q} | \
      perl -lape 's/graph\[/graph[rankdir=LR,/g' | \
      dot -Tsvg > {output:q}
    '''

rule rulegraph_svg:
    input: 'Snakefile', 'rulegraph.Snakefile'
    output: 'rulegraph-{target}.svg'
    shell: '''
    set -o pipefail;
    snakemake --nolock -f --rulegraph {wildcards.target:q} | \
      perl -lape 's/graph\[/graph[rankdir=TB,/g' | \
      dot -Tsvg > {output:q}
    '''
