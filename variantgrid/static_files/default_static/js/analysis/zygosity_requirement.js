/* The zygosity table in the Trio / Duo / Quad node editors: what one member's options read as under a
   ZygosityRequirement (patients.enums.ZygosityRequirement codes). REF stands for REF or no record unless
   a REF call is required, so that case says which it is - the same rule as _zygosity_options in Python.
   Node editors load by ajax, so this file is evaluated once per editor opened: it declares only functions,
   as a top-level const can't be declared twice */
function zygosityRequirementOptions(options, requirement) {
    const ALLOW_NO_CALL = 'N';
    const REF_CALL_REQUIRED = 'R';

    if (requirement === ALLOW_NO_CALL) {
        return options + ', ?';
    }
    if (requirement === REF_CALL_REQUIRED) {
        return options.replace(/\bREF\b/, 'REF (called)');
    }
    return options;
}
