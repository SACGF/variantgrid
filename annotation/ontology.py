class OntologyMixin:

    @property
    def id_str(self) -> str:
        return self.int_as_id(self.pk)

    @property
    def url(self) -> str:
        # FIXME redundant to dbregex but don't want to reference that from this file
        raise NotImplementedError(f"Ontology {self} has not defined a URL method")

    @classmethod
    def int_as_id(cls, ontology_index: int) -> str:
        num_part = str(ontology_index).rjust(cls.expected_length(), '0')
        return f"{cls.PREFIX}{num_part}"

    @classmethod
    def id_as_int(cls, ontology_id_str: str) -> str:
        parts = ontology_id_str.split(":")
        if len(parts) != 2:
            raise ValueError(f"Illegal Ontology String ({ontology_id_str})")
        # somewhat redundant that the prefixes include the ":" character
        if parts[0] != cls.PREFIX[:-1]:
            raise ValueError(f"This doesn't appear to be the right id type - expected {cls.PREFIX[:-1]} got {parts[0]}")

        return int(parts[1])

    @classmethod
    def expected_length(cls) -> int:
        raise NotImplementedError(f"Ontology {cls} has not defined expected length")

    @property
    def padded_id(self) -> str:
        return str(self.pk).rjust(self.expected_length(), '0')
