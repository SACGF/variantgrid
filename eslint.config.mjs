import globals from "globals";

export default [
    {
        // global ignores (must be a standalone object to apply to all configs)
        ignores: [
            "**/*.min.js",
            "**/js/lib/**", // vendored third-party libraries
            "variantgrid/sitestatic/**", // collectstatic output
            "node_modules/**",
        ],
    },
    {
        files: ["variantgrid/static_files/**/*.js"],
        languageOptions: {
            ecmaVersion: 2022,
            sourceType: "script", // classic scripts loaded via <script src>, not ES modules
            globals: {
                ...globals.browser,
                ...globals.jquery,
                // vendored libraries loaded via <script> tags
                jsPlumb: "readonly",
                moment: "readonly",
                Plotly: "readonly",
                Cookies: "readonly",
                Split: "readonly",
                Rollbar: "readonly",
                chainedfk: "readonly", // django-smart-selects
                chainedm2m: "readonly", // django-smart-selects
                Urls: "readonly", // django-js-reverse
                VGLoaders: "readonly", // js/lib/loaders/vg-loaders.js
                // cross-file project globals (defined at top level of first-party files)
                EKey: "readonly", // vc_keys.js
                EKeys: "readonly", // vc_keys.js
                SpecialEKeys: "readonly", // vc_keys.js
                VcSettings: "readonly", // vc_settings.js
                CitationsManager: "readonly", // citations.js
                createModalShell: "readonly", // global.js
                showReloadPageErrorDialog: "readonly", // global.js
                limitLength: "readonly", // global.js
                debounce: "readonly", // global.js
                escapeHtml: "readonly", // global.js
                svgSelect: "readonly", // global.js
                toFixedString: "readonly", // scientific_number_widget.js
                toPercent: "readonly", // scientific_number_widget.js
                loadNodeData: "readonly", // analysis.js
                getGridAndEditorWindow: "readonly", // analysis.js
                lockNodeField: "readonly", // analysis.js
                replaceEditorWindow: "readonly", // analysis.js
                loadGridAndEditorForNode: "readonly", // analysis.js
                analysisVariable: "readonly", // analysis.js
                addAnalysisVariableButton: "readonly", // analysis.js
                _getAnalysisWindow: "readonly", // analysis.js
                createJSEvent: "readonly", // analysis.js
                addVariantTag: "readonly", // analysis.js
                removeVariantTag: "readonly", // analysis.js
                hideLoadingOverlay: "readonly", // analysis.js
                showBottomPaneTab: "readonly", // analysis.js
                BOTTOM_PANE_GRID: "readonly", // analysis.js
                openVariantDetailsTab: "readonly", // analysis.js
                hideGridLoadingOverlay: "readonly", // analysis.js
                getAnalysisWindow: "readonly", // grid.js
                createGridLink: "readonly", // grid.js
                inAnalysis: "readonly", // grid.js
                getVariantTagHtml: "readonly", // grid.js
                sortVariantTags: "readonly", // grid.js
                variantTaggingPillOptions: "readonly", // grid.js
                createIgvUrl: "readonly", // igv.js (typeof-guarded use in vc_links.js)
                createIgvLink: "readonly", // igv.js
                renderIgvLocusLinks: "readonly", // igv.js
                deleteNodesFromDOM: "readonly", // analysis_nodes.js
                unselectActive: "readonly", // analysis_nodes.js
                addNodesToDOM: "readonly", // analysis_nodes.js
                attatchAnalysisNodeConnections: "readonly", // analysis_nodes.js
                checkAndMarkDirtyNodes: "readonly", // analysis_nodes.js
                getNode: "readonly", // analysis_nodes.js
                loggedOutHandler: "readonly", // analysis_nodes.js
                createSampleNode: "readonly", // samplenode.js
                poll_cached_generated_file: "readonly", // cached_generated_files.js
                venn2: "readonly", // venn_intersect.js
                venn_select: "readonly", // venn_intersect.js
                renderNodeIcon: "readonly", // node_type_select.js
                buildTipBox: "readonly", // tips.js
                addGeneEvidence: "readonly", // panel_app.js
                getPanelAppGeneEvidenceDiv: "readonly", // panel_app.js
                loadPatientPhenotype: "readonly", // patient_phenotypes.js
                patientPhenotypeHtml: "readonly", // patient_phenotypes.js
                getOntologyTermLinks: "readonly", // phenotype.js
                expandCollapsedOntologyTerm: "readonly", // phenotype.js
                defaultLayout: "readonly", // plotly_helper.js
                plotHBarArrays: "readonly", // plotly_helper.js
                SampleSelectionActions: "readonly", // sample_selection_actions.js
                VariantGridFilterBuilder: "readonly", // variantgrid_filter_builder.js
                VariantGridFormat: "readonly", // variantgrid_formats.js
                getValue: "readonly", // global.js
                removeItemFromArray: "readonly", // global.js
                checkLoggedIn: "readonly", // global.js
                convertTimestamp: "readonly", // global.js
                JS_DATE_FORMAT_SECONDS: "readonly", // global.js
                JS_DATE_FORMAT_SCIENTIFIC: "readonly", // global.js
                loadAjaxBlock: "readonly", // global.js
                formatJson: "readonly", // global.js
                createTimestampDom: "readonly", // global.js
                setupModalAnimationForWebTesting: "readonly", // global.js
                EncodeQueryData: "readonly", // global.js
                deleteItemClickHandler: "readonly", // global.js
                createModal: "readonly", // global.js
                highlightTextAsDom: "readonly", // global.js
                limitLengthSpan: "readonly", // global.js
                FloatingPanel: "readonly", // global.js
                closeSelect2Dropdowns: "readonly", // global.js
                clearAutocompleteChoice: "readonly", // global.js
                dynamicSort: "readonly", // global.js
                createMessage: "readonly", // global.js
                DataTableDefinition: "readonly", // datatable_definition.js
                VCTable: "readonly", // vc_form.js
                Flags: "readonly", // flags.js
                VCLinks: "readonly", // vc_links.js
                // cross-execution state in fragment-loaded files: these scripts are
                // re-executed on AJAX load, so the state must stay a plain assignment
                // (a top-level let/const would throw on re-execution)
                RAISED_GET_ANALYSIS_WINDOW_JS_ERROR: "writable", // grid.js
                VENN_TOGGLE_WIDGET_CLASS: "writable", // venn_intersect.js
                venn_id: "writable", // venn_intersect.js
                freq: "writable", // cached_generated_files.js
                // globals injected by Django templates (inline <script> blocks)
                ANALYSIS_ID: "readonly", // analysis.html
                ANALYSIS_SETTINGS: "readonly", // analysis_settings_node_counts_tab.html
                ANALYSIS_HORIZONTAL_MODE: "readonly", // analysis.html
                ANALYSIS_LOADING_ANIMATIONS: "readonly", // analysis.html
                ANALYSIS_TIPS: "readonly", // analysis.html
                ANALYSIS_TAGS_NODE_ID: "readonly", // analysis.html
                NODE_HELP: "readonly", // analysis.html
                messagePoller: "readonly", // analysis.html
                analysisNodeVariables: "readonly", // analysis.html
                saveSettingsOnResize: "writable", // analysis.html, assigned in analysis.js
                secondWindow: "writable", // analysis.html, assigned in analysis.js
                panelResizeTimeout: "writable", // analysis.html, assigned in analysis.js
                panelResizeUpdateDelay: "readonly", // analysis.html
                variantTags: "readonly", // analysis.html
                loadInitialGridEditor: "readonly", // analysis.html
                registerComponent: "readonly", // analysis_editor_and_grid.html
                EDITOR: "readonly", // analysis_editor_and_grid.html
                reloadNodes: "readonly", // analysis_settings.html
                variantTagStaleDays: "writable", // analysis.html, reassigned in analysis_nodes.js
                nodeProbandSampleId: "readonly", // node_data_grid.html, sample_variants_tab.html
                nodeProbandPatientId: "readonly", // node_data_grid.html
            },
        },
        rules: {
            "no-var": "error",
            "prefer-const": ["error", {"destructuring": "all"}],
            "no-undef": "warn",
            // builtinGlobals false: files that define a shared global (e.g. const EKeys
            // in vc_keys.js) would otherwise be flagged for redeclaring the config entry
            "no-redeclare": ["error", {"builtinGlobals": false}],
            "semi": ["error", "always"],
        },
    },
];
