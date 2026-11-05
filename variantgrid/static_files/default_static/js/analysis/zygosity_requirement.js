/* The zygosity table in the Trio / Duo / Quad node editors: what one member's options read as under a
   ZygosityRequirement (patients.enums.ZygosityRequirement codes). REF stands for REF or no record unless
   a REF call is required, so that case says which it is - the same rule as _zygosity_options in Python */
const ZYGOSITY_REQUIREMENT_ALLOW_NO_CALL = 'N';
const ZYGOSITY_REQUIREMENT_REF_CALL_REQUIRED = 'R';

function zygosityRequirementOptions(options, requirement) {
    if (requirement === ZYGOSITY_REQUIREMENT_ALLOW_NO_CALL) {
        return options + ', ?';
    }
    if (requirement === ZYGOSITY_REQUIREMENT_REF_CALL_REQUIRED) {
        return options.replace(/\bREF\b/, 'REF (called)');
    }
    return options;
}
