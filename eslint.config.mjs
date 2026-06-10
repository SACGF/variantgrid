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
                getAnalysisWindow: "readonly", // grid.js
                deleteNodesFromDOM: "readonly", // analysis_nodes.js
                unselectActive: "readonly", // analysis_nodes.js
                createSampleNode: "readonly", // samplenode.js
                // globals injected by Django templates (inline <script> blocks)
                ANALYSIS_ID: "readonly", // analysis.html
                ANALYSIS_SETTINGS: "readonly", // analysis_settings_node_counts_tab.html
                NODE_HELP: "readonly", // analysis.html
                messagePoller: "readonly", // analysis.html
                analysisNodeVariables: "readonly", // analysis.html
                saveSettingsOnResize: "writable", // analysis.html, assigned in analysis.js
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
