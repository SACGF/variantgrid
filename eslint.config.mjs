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
            },
        },
        rules: {
            "no-var": "error",
            "prefer-const": ["error", {"destructuring": "all"}],
            "no-undef": "warn",
            "no-redeclare": "error",
        },
    },
];
