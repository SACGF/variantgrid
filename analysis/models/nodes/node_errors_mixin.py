""" NodeErrorsMixin - what is wrong with an AnalysisNode: analysis, parent, internal (load) and
    configuration errors as (NodeErrorSource, message) pairs, field errors the user may waive, warnings,
    and the node status / exception an error source maps to. Entry points: get_errors, get_warnings,
    get_field_errors, get_status_from_errors, throw_errors_exception. """
from analysis.exceptions import (
    NodeConfigurationException,
    NodeParentErrorsException,
    NonFatalNodeError,
)
from analysis.models.enums import AnalysisTemplateType, NodeErrorSource, NodeStatus
from snpdb.models import GenomeBuild


class NodeErrorsMixin:
    """ Part of AnalysisNode. Hooks a node subclass overrides:
        - _get_configuration_errors: messages for a node that can't run as configured
        - _get_field_errors: errors keyed by the field at fault, shown against it in the editor;
          get_ignorable_error_fields names those ignore_field_errors may waive (then reported as warnings)
        - get_warnings: shown on the node without stopping it - extend via super() """

    def _get_genome_build_errors(self, field_name, field_genome_build: GenomeBuild) -> list:
        """ Used to quickly add errors about genome build mismatches
            This only happens in templates (ran template on sample with different build than hardcoded data)
            In normal analyses, autocomplete restrictions should not allow you to configure data from other builds """
        errors = []
        if field_genome_build != self.analysis.genome_build:
            msg = f"{field_name} genome build: {field_genome_build} different from analysis build: {self.analysis.genome_build}"
            errors.append(msg)
        return errors

    def _get_configuration_errors(self) -> list:
        return []

    def get_parents_and_errors(self):
        """ Returns error array, includes any min/max parent error and node config error """
        if self.has_input():
            return self.get_parent_subclasses_and_errors()
        return [], []

    def _get_analysis_errors(self) -> list[str]:
        if self._cached_analysis_errors is None:
            self._cached_analysis_errors = self.analysis.get_errors()
        return self._cached_analysis_errors

    def get_analysis_errors(self, flat=False):
        """ Analysis level errors (eg no annotation version) - these can't be fixed by editing the node """
        errors = [(NodeErrorSource.ANALYSIS, e) for e in self._get_analysis_errors()]
        if flat:
            errors = NodeErrorsMixin.flatten_errors(errors)
        return errors

    def get_ignorable_error_fields(self) -> list[str]:
        """ Fields whose _get_field_errors() the user may waive with ignore_field_errors """
        return []

    def _get_field_errors(self) -> dict[str, list[str]]:
        """ Configuration errors keyed by the node field at fault, so the editor shows them against
            that field. Templates are exempt - they're configured however the author likes """
        return {}

    def _get_ignored_error_fields(self) -> set[str]:
        if self.ignore_field_errors:
            return set(self.get_ignorable_error_fields())
        return set()

    def get_field_errors(self) -> dict[str, list[str]]:
        """ The field errors that stop the node running """
        if self.analysis.template_type == AnalysisTemplateType.TEMPLATE:
            return {}
        ignored = self._get_ignored_error_fields()
        return {field: errors for field, errors in self._get_field_errors().items()
                if errors and field not in ignored}

    def get_ignored_field_errors(self) -> list[str]:
        if self.analysis.template_type == AnalysisTemplateType.TEMPLATE:
            return []
        ignored = self._get_ignored_error_fields()
        return [e for field, errors in self._get_field_errors().items() if field in ignored for e in errors]

    def get_warnings(self) -> list[str]:
        return list(self.get_ignored_field_errors())

    def _get_non_field_errors(self, include_parent_errors=True) -> list[tuple[NodeErrorSource, str]]:
        errors = self.get_analysis_errors()
        _, parent_errors = self.get_parents_and_errors()
        if include_parent_errors:
            errors.extend(parent_errors)
        if self.errors:
            errors.append((NodeErrorSource.INTERNAL_ERROR, self.errors))
        errors.extend((NodeErrorSource.CONFIGURATION, ce) for ce in self._get_configuration_errors())
        return errors

    def get_errors(self, include_parent_errors=True, flat=False):
        """ returns a tuple of (NodeError, str) unless flat=True where it's only string """
        errors = self._get_non_field_errors(include_parent_errors=include_parent_errors)
        errors.extend((NodeErrorSource.CONFIGURATION, e) for field_errors in self.get_field_errors().values()
                      for e in field_errors)
        if flat:
            errors = NodeErrorsMixin.flatten_errors(errors)
        return errors

    def can_ignore_errors(self) -> bool:
        """ Everything wrong with the node is a field error ignore_field_errors would waive - what the
            Template tab offers when it reveals a branch the template run hid """
        field_errors = self.get_field_errors()
        if not field_errors or self._get_non_field_errors():
            return False
        return set(field_errors) <= set(self.get_ignorable_error_fields())

    @staticmethod
    def flatten_errors(errors):
        return [f"{NodeErrorSource(nes).label}: {error}" for nes, error in errors]

    def get_error_summaries(self) -> list[str]:
        """ get_errors(flat=True) for the tooltip on the node card's error cross - a load's traceback
            isn't for users, so it's only named as an internal error """
        return [NodeErrorSource.INTERNAL_ERROR.label if source == NodeErrorSource.INTERNAL_ERROR
                else f"{NodeErrorSource(source).label}: {error}"
                for source, error in self.get_errors()]

    @staticmethod
    def get_status_from_errors(errors):
        ERROR_STATUS = {
            NodeErrorSource.INTERNAL_ERROR: NodeStatus.ERROR,
            NodeErrorSource.ANALYSIS: NodeStatus.ERROR_WITH_PARENT,
            NodeErrorSource.PARENT: NodeStatus.ERROR_WITH_PARENT,
            NodeErrorSource.CONFIGURATION: NodeStatus.ERROR_CONFIGURATION,
        }
        if not errors:
            raise ValueError("Passed in empty errors!")

        error_sources = {s for s, _ in errors}
        for source, status in ERROR_STATUS.items():
            if source in error_sources:
                return status
        raise ValueError("No error source found")

    @staticmethod
    def throw_errors_exception(errors):
        ERROR_EXCEPTIONS = {
            NodeErrorSource.INTERNAL_ERROR: ValueError,
            NodeErrorSource.ANALYSIS: NonFatalNodeError,
            NodeErrorSource.PARENT: NodeParentErrorsException,
            NodeErrorSource.CONFIGURATION: NodeConfigurationException,
        }
        if not errors:
            raise ValueError("Passed in empty errors!")
        error_sources = {s for s, _ in errors}
        for source, exception_klass in ERROR_EXCEPTIONS.items():
            if source in error_sources:
                raise exception_klass()
        raise ValueError("No error source found")
