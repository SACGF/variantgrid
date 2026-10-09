// @ts-check
// snpdb/templates/snpdb/patients/cohort_gene_counts.html
$(document).ready(function() {
    const data = readJsonData("cohort-gene-counts-data");
    const gcContainer = $("#cohort-gene-counts-graph-container");

    function load_graph(url) {
        gcContainer.html('<i class="fa fa-spinner"></i>');
        gcContainer.load(url);
    }

    let geneListId = data.gene_list_id;

    function load_gene_list() {
        const geneCountType = $("#id_gene_count_type").val();
        if (geneCountType && geneListId) {
            load_graph(Urls.cohort_gene_counts_matrix(data.cohort_id, geneCountType, geneListId));
        } else {
            gcContainer.empty();
        }
    }

    $('#id_gene_count_type').change(load_gene_list);

    $('#id_gene_list').change(function() {
        geneListId = $(this).val();
        load_gene_list();
    });

    load_gene_list();

});
