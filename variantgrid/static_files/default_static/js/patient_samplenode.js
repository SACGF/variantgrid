function createExampleSamplenode(patientDict) {
	const container = $("#sample-node-container");
	container.empty();
	const nodeData = {name : patientDict['name'], attributes: {'class' : 'window'}};
	const sample_node = createSampleNode(nodeData);
	sample_node.each(function() { this.updateState({patient: patientDict}); });
	container.append(sample_node);
}


function getName() {
    const last_name = $("#id_last_name").val();
    const first_name = $("#id_first_name").val();
    let full_name = '';
    if (first_name) {
        full_name += first_name;         
    }

    if (last_name) {
        if (full_name) {
            full_name += '\n';                
        }
        full_name += last_name;
    }
    return full_name;
}


function getPatientDict(form) {
	const name = getName();
	const sex = $('#id_sex', form).val();
	const deceased = !!$("#id_date_of_death").val();
	return {name: name, sex : sex, deceased: deceased};
}

function formChanged() {
	const patientDict = getPatientDict(this);
	createExampleSamplenode(patientDict);
}

function setupPatientSamplenode(initialPatientDict) {
	const form = $("form#patient-form");

	// Create a clear / init / reset function
	form.each(function() {
		const resetPatientSampleNode = function() {
			createExampleSamplenode(initialPatientDict);
		};
		this.resetPatientSampleNode = resetPatientSampleNode;
		resetPatientSampleNode(); // set initial
	});
	form.change(formChanged);
	$("input[type=submit]", form).button();
	$("button#reset-form").button().click(function() {
		form.each(function() {
			this.reset();
			this.resetPatientSampleNode();
		});
	});
}

