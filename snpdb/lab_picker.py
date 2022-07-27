import itertools
from dataclasses import dataclass
from typing import Optional, Union, List, Set, Iterator, Iterable

from django.contrib.auth.models import User
from django.db.models import QuerySet
from django.http import HttpRequest, HttpResponseRedirect
from django.shortcuts import redirect
from django.urls import reverse
from lazy import lazy
from snpdb.models import Lab, Organization, GenomeBuild, UserSettings


@dataclass
class OrgLabGroup:
    org: Organization
    labs: List[Lab]
    external: bool = False

    @property
    def is_multi_lab(self):
        return len(self.labs) > 1


@dataclass(frozen=True)
class LabSelection:
    """
    Used to simply LabPickerData, tracks which labs the view has access to, and which ones are selected
    """
    all_labs: List[Lab]
    selected_labs: Set[Lab]
    selected_org: Optional[Organization]
    selected_all: bool

    @staticmethod
    def from_user(user: User, selection: Union[str, int]) -> 'LabSelection':
        """
        Provide a user, a lab selection string. Method ensures use has access to the selection
        :param user: The user
        :param selection: primary key of a lab, or "org-{org_id}" or 0 or blank to indicate all labs the user can access
        :return: A LabSelection
        """
        all_labs = Lab.valid_labs_qs(user=user, admin_check=True).exclude(organization__active=False).order_by(
            'organization__name', 'name')
        seleced_labs: QuerySet[Lab]
        selected_org: Optional[Organization] = None
        selected_all = False

        if isinstance(selection, str) and selection.isnumeric():
            try:
                # check to see if selection is an
                selection = int(selection)
            except ValueError:
                pass

        if isinstance(selection, int) and not selection == 0:
            # legacy method of just passing LabID around
            selected_labs = all_labs.filter(pk=selection)
        elif isinstance(selection, str) and selection.startswith("org-"):
            org_pk = int(selection[4:])
            selected_org = Organization.objects.filter(pk=org_pk).get()
            selected_labs = all_labs.filter(organization=selected_org)
        elif isinstance(selection, Lab):
            # make sure the user has access to the passed in lab
            selected_labs = all_labs.filter(pk=selection.pk)
        else:
            selected_labs = all_labs.filter(external=False)
            selected_all = True

        if not selected_labs.exists():
            raise ValueError(f"You do not have access to lab selection : {selection}")

        return LabSelection(
            all_labs=list(all_labs),
            selected_labs=set(selected_labs),
            selected_org=selected_org,
            selected_all=selected_all
        )

    @staticmethod
    def single_lab(lab: Lab) -> 'LabSelection':
        """
        Use just when the user is not relevant
        :param lab: The lab to be selected, no checks performed
        """
        return LabSelection(
            all_labs=[lab],
            selected_labs={lab, },
            selected_org=None,
            selected_all=False
        )

    @staticmethod
    def multi_lab(labs: Iterable[Lab]) -> 'LabSelection':
        """
       Use just when the user is not relevant
       :param labs: The lab to be selected, no checks performed
       """
        labs = list(labs)
        return LabSelection(
            all_labs=labs,
            selected_labs=set(labs),
            selected_org=None,
            selected_all=False
        )

    @property
    def selection(self) -> Union[str, int]:
        """
        Generates a value that can be used to generate another LabSelection for the same labs.
        """
        if self.selected_all:
            return 0
        if self.selected_org:
            return f"org-{self.selected_org.pk}"
        return ",".join([f"{lab.pk}" for lab in self.selected_labs])


@dataclass(frozen=True)
class LabPickerData:
    """
    TODO: This class has taken over from UserPerspective, maybe it should be renamed?
    Useful for managing: pages with lab pickers
    Methods which can operate on one or more labs: e.g. allele graph
    Helper methods for getting genome build
    Can be user focused or lab focused (e.g. for lab wide notifications)
    """

    lab_selection: LabSelection  # what labs are available to the user, and what ones have been selected by the user
    user: Optional[User] = None  # the user who we're seeing the POV of
    view_name: Optional[str] = None  # useful for when the page needs to be reloaded if lab changes
    multi_select: bool = True  # can we select multiple labs on this page

    @property
    def selection(self):
        return self.lab_selection.selection

    @property
    def all_labs(self) -> List[Lab]:
        return self.lab_selection.all_labs

    @property
    def selected_labs(self) -> Set[Lab]:
        return self.lab_selection.selected_labs

    @property
    def selected_org(self) -> Optional[Organization]:
        return self.lab_selection.selected_org

    @property
    def internal_labs(self):
        return [lab for lab in self.all_labs if not lab.external]

    @property
    def external_labs(self):
        return [lab for lab in self.all_labs if lab.external]

    @staticmethod
    def from_request(
        request: HttpRequest,
        selection: Optional[Union[str, int]] = None,
        view_name: Optional[str] = None,
        multi_select: bool = True
    ) -> 'LabPickerData':
        """
        Just calls for_user with request.user, but could do more in future
        """
        return LabPickerData.for_user(
            user=request.user,
            selection=selection,
            view_name=view_name,
            multi_select=multi_select
        )

    @staticmethod
    def for_user(
            user: User,
            selection: Optional[Union[str, int]] = None,
            view_name: Optional[str] = None,
            multi_select: bool = True
    ) -> 'LabPickerData':
        return LabPickerData(
            lab_selection=LabSelection.from_user(user=user, selection=selection),
            user=user,
            view_name=view_name,
            multi_select=multi_select
        )

    @staticmethod
    def for_lab(lab: Lab, user: Optional[User] = None):
        return LabPickerData(
            lab_selection=LabSelection.single_lab(lab=lab),
            user=user
        )

    @staticmethod
    def for_labs(labs: Iterable[Lab], user: Optional[User] = None):
        return LabPickerData(
            lab_selection=LabSelection.multi_lab(labs=labs),
            user=user
        )

    @lazy
    def genome_build(self) -> GenomeBuild:
        return UserSettings.get_for(lab=self.selected_lab, organization=self.selected_org).default_genome_build or GenomeBuild.grch38()

    @property
    def your_labs(self):
        # provided as compatibility for UserPerspective
        return set(self.all_labs)

    @property
    def labs_if_not_admin(self) -> Set[Lab]:
        # Used to group together as "internal" labs
        if self.user and self.user.is_superuser:
            return set()
        else:
            return set(self.all_labs)

    @property
    def selected_lab(self) -> Optional[Lab]:
        if len(self.selected_labs) == 1:
            return list(self.selected_labs)[0]
        return None

    @property
    def multi_labs_selected(self) -> bool:
        return len(self.selected_labs) > 1

    @property
    def has_multi_labs(self) -> bool:
        return len(self.internal_labs) + len(self.external_labs) > 1

    @property
    def has_multi_orgs(self) -> bool:
        if self.all_labs:
            first_org = self.all_labs[1].organization
            for lab in self.all_labs[1:]:
                if lab.organization != first_org:
                    return True
        return False

    @property
    def has_external(self) -> bool:
        return bool(self.external_labs)

    def is_admin_mode(self) -> bool:
        # see if user is an admin looking at all labs
        return self.user and self.user.is_superuser and (not self.selected_lab) and (not self.selected_org)

    def org_lab_groups(self) -> Iterator[OrgLabGroup]:
        for org, labs in itertools.groupby(self.internal_labs, lambda l: l.organization):
            yield OrgLabGroup(org, list(labs))
        for org, labs in itertools.groupby(self.external_labs, lambda l: l.organization):
            yield OrgLabGroup(org, list(labs), external=True)

    def check_redirect(self) -> Optional[HttpResponseRedirect]:
        """
        LabPicker pages redirect to a specific URL if only one lab is available, e.g. if
        /classifications/foo (shows all labs, but the user only has access to lab 3) calling this will generate
        /classifications/foo/3
        Also (and more importantly) for pages that don't support multi_select
        :return:
        """
        if len(self.all_labs) == 1 and self.selection is None:
            return redirect(reverse(self.view_name, kwargs={'lab_id': self.selected_lab.pk}))
        if not self.multi_select and not self.selected_lab:  # need only 1 lab selected in this case
            default_lab = UserSettings.get_for_user(self.user).default_lab_safe()
            return redirect(reverse(self.view_name, kwargs={'lab_id': default_lab.pk}))

    @property
    def lab_ids(self) -> Set[int]:
        """
        :return: A set of selected lab.pks
        """
        return {lab.pk for lab in self.selected_labs}
