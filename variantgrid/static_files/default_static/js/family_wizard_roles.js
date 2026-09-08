/* Trio/Quad wizards - keeps the pedigree figure in step with the roles and affected ticks the user
   picks. Wants one select per sample inside '.sample-select' and, inside '.sample-affected', an
   '#id_sample_<n>_affected' checkbox plus the '.proband-affected' note that stands in for it on the
   proband's row - all in the same order as sampleSexes.

   roleClasses maps a role value to the class that fills that symbol in (see --pedigree-*-fill in
   uicore/tags/svg_icon_sprite.html). autoAssignRoles maps a sex to the role a sample of that sex
   takes once the proband is chosen - pass it where sex, plus what is left over, decide the rest. */
const FAMILY_PROBAND = 'P';

function setupFamilyRoles(sampleSexes, roleClasses, autoAssignRoles) {
    const roleSelects = $("select", ".sample-select");

    function affectedCheckbox(i) {
        return $("#id_sample_" + (i + 1) + "_affected");
    }

    function setRole(i, role) {
        const select = roleSelects.eq(i);
        if (select.find("option[value='" + role + "']").length) {
            select.val(role);
        }
    }

    function autoAssignFromProband(probandIndex) {
        let samplesLeft = [];
        roleSelects.each(function(i) {
            if (i !== probandIndex) {
                samplesLeft.push(i);
            }
        });
        let rolesLeft = Object.values(autoAssignRoles);

        // Sex settles a sample's role, unless two of them claim the same one
        for (const [sex, role] of Object.entries(autoAssignRoles)) {
            const claimants = samplesLeft.filter(i => sampleSexes[i]["sex"] === sex);
            if (claimants.length === 1) {
                setRole(claimants[0], role);
                samplesLeft = samplesLeft.filter(i => i !== claimants[0]);
                rolesLeft = rolesLeft.filter(r => r !== role);
            }
        }

        // One sample and one role left over take each other
        if (samplesLeft.length === 1 && rolesLeft.length === 1) {
            setRole(samplesLeft[0], rolesLeft[0]);
        }
    }

    function refresh() {
        const filled = {};
        roleSelects.each(function(i) {
            const role = $(this).val();
            const checkbox = affectedCheckbox(i);
            const isProband = role === FAMILY_PROBAND;
            // The proband is affected by definition, so say so rather than offering a switch
            checkbox.prop("disabled", isProband).prop("checked", isProband || checkbox.is(":checked"));
            checkbox.closest(".custom-switch").toggle(!isProband);
            checkbox.closest(".sample-affected").find(".proband-affected").toggle(isProband);
            const cssClass = roleClasses[role];
            if (cssClass) {
                filled[cssClass] = checkbox.is(":checked");
            }
        });
        for (const cssClass of Object.values(roleClasses)) {
            $(".pedigree-figure").toggleClass(cssClass, Boolean(filled[cssClass]));
        }
    }

    roleSelects.change(function() {
        if (autoAssignRoles && $(this).val() === FAMILY_PROBAND) {
            autoAssignFromProband(roleSelects.index(this));
        }
        refresh();
    });
    $("input", ".sample-affected").change(refresh);
    refresh();
}
