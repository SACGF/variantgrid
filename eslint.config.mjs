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
                _: "readonly",
                d3: "readonly",
                Plotly: "readonly",
                Cookies: "readonly",
                Split: "readonly",
                Rollbar: "readonly",
                chainedfk: "readonly", // django-smart-selects
                chainedm2m: "readonly", // django-smart-selects
                Urls: "readonly", // django-js-reverse
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
                getAnalysisWindow: "readonly", // grid.js
                createIgvUrl: "readonly", // grid.js (typeof-guarded use in vc_links.js)
                deleteNodesFromDOM: "readonly", // analysis_nodes.js
                unselectActive: "readonly", // analysis_nodes.js
                addNodesToDOM: "readonly", // analysis_nodes.js
                attatchAnalysisNodeConnections: "readonly", // analysis_nodes.js
                checkAndMarkDirtyNodes: "readonly", // analysis_nodes.js
                getNode: "readonly", // analysis_nodes.js
                loggedOutHandler: "readonly", // analysis_nodes.js
                createSampleNode: "readonly", // samplenode.js
                createTrioNode: "readonly", // pedigree_node.js
                createQuadNode: "readonly", // pedigree_node.js
                venn2: "readonly", // venn_intersect.js
                venn_select: "readonly", // venn_intersect.js
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
                DataTableDefinition: "readonly", // datatable_definition.js
                VCTable: "readonly", // vc_form.js
                Flags: "readonly", // flags.js
                VCLinks: "readonly", // vc_links.js
                // cross-execution state in fragment-loaded files: these scripts are
                // re-executed on AJAX load, so the state must stay a plain assignment
                // (a top-level let/const would throw on re-execution)
                seen_igv_error: "writable", // grid.js
                RAISED_GET_ANALYSIS_WINDOW_JS_ERROR: "writable", // grid.js
                VENN_TOGGLE_WIDGET_CLASS: "writable", // venn_intersect.js, read by analysis_nodes.js
                venn_id: "writable", // venn_intersect.js
                freq: "writable", // cached_generated_files.js
                // globals injected by Django templates (inline <script> blocks)
                ANALYSIS_ID: "readonly", // analysis.html
                ANALYSIS_SETTINGS: "readonly", // analysis_settings_node_counts_tab.html
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
