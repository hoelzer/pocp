process plot {
    label 'python'
    publishDir "${params.output}", mode: 'copy', pattern: "pocp-heatmap.{svg,pdf}"

    input:
      path(pocp_matrix)

    output:
	    path('pocp-heatmap.svg')
	    path('pocp-heatmap.pdf')

    script:
    """
    plot-heatmap.py --matrix ${pocp_matrix} --width ${params.width} --height ${params.height}
    """

    stub:
    """
    touch pocp-heatmap.svg pocp-heatmap.pdf
    """
}
