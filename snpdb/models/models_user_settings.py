import dataclasses
from dataclasses import dataclass

from avatar.templatetags.avatar_tags import avatar_url
from collections import OrderedDict, defaultdict

from django.contrib.auth.models import User, Group
from django.core.exceptions import ValidationError
from django.db import models
from django.db.models.deletion import SET_NULL, CASCADE
from django_extensions.db.models import TimeStampedModel
from typing import Optional, List, Tuple, Dict, Set

from lazy import lazy
from model_utils.managers import InheritanceManager

from library.django_utils import thread_safe_unique_together_get_or_create
from library.utils import rgb_invert, string_deterministic_hash
from snpdb.models.models import Tag, Lab, Organization
from snpdb.models.models_columns import CustomColumnsCollection, CustomColumn
from snpdb.models.models_genome import GenomeBuild
from snpdb.models.models_enums import BuiltInFilters


def get_igv_data(user, genome_build: GenomeBuild = None):
    user_settings = UserSettings.get_for_user(user)
    replace_dict = UserDataPrefix.get_replace_dict(user)
    igv_data = {'base_url': "http://localhost:%d" % user_settings.igv_port,
                'replace_dict': replace_dict}
    if genome_build:
        igv_data['genome'] = genome_build.igv_genome
    return igv_data


class UserDataPrefix(models.Model):
    """ This is used to alter user paths (eg from /data/x on the server to a mapped network drive etc) """
    user = models.ForeignKey(User, on_delete=CASCADE)
    prefix = models.TextField()
    replacement = models.TextField()

    @staticmethod
    def get_replace_dict(user):
        replace_dict = {}
        for prefix, replacement in user.userdataprefix_set.values_list("prefix", "replacement"):
            replace_dict[prefix] = replacement
        return replace_dict


class UserTagColors(models.Model):
    user = models.ForeignKey(User, null=True, on_delete=CASCADE)  # Null is "default"
    tag = models.ForeignKey(Tag, on_delete=CASCADE)
    rgb = models.CharField(max_length=7)  # '#rrggbb'

    @staticmethod
    def get_all_tags_and_user_colors(user):
        """ Returns Hash of { tag_name : color }
            with color being None if not set for user """
        user_colors_by_tag = {}
        for tag_id, rgb in UserTagColors.objects.filter(user=user).values_list('tag', 'rgb'):
            user_colors_by_tag[tag_id] = {"background-color": rgb,
                                          "color": rgb_invert(rgb)}

        user_tag_colors = OrderedDict()
        for tag in Tag.objects.all():
            user_tag_colors[tag] = user_colors_by_tag.get(tag.id)
        return user_tag_colors

    @staticmethod
    def get_tag_styles_and_colors(user):
        user_tag_styles = []
        user_tag_colors = {}

        for tag, style in UserTagColors.get_all_tags_and_user_colors(user).items():
            user_tag_styles.append((tag.id, style))
            rgb = None
            if style:
                rgb = style.get("background-color")
            user_tag_colors[tag.id] = rgb
        return user_tag_styles, user_tag_colors


class UserPageAck(TimeStampedModel):
    """ Used for eg I've read the notices on this page """
    user = models.ForeignKey(User, on_delete=CASCADE)
    page_id = models.TextField()

    class Meta:
        unique_together = ("user", "page_id")


class UserGridConfig(models.Model):
    user = models.ForeignKey(User, on_delete=CASCADE)
    grid_name = models.TextField()  # JQGrid caption
    rows = models.IntegerField(default=10)
    show_group_data = models.BooleanField(default=True)
    show_incomplete_data = models.BooleanField(default=False)
    show_hidden_data = models.BooleanField(default=False)
    filter_name = models.TextField(null=True)  # Can use this to make arbitrary filters

    class Meta:
        unique_together = ("user", "grid_name")

    @staticmethod
    def get(user: User, grid_name: str) -> 'UserGridConfig':
        return thread_safe_unique_together_get_or_create(UserGridConfig, user=user, grid_name=grid_name)[0]

    def __str__(self):
        details = []
        if self.show_group_data:
            details.append("group")
        if self.show_incomplete_data:
            details.append("incomplete")
        if self.show_hidden_data:
            details.append("hidden")

        return f"{self.user}/{self.grid_name}: {','.join(details)}"


class SettingsOverride(models.Model):
    """ We want UserSettings to cascade via Org -> Lab -> User
        Where the later levels override if they are non-null """
    objects = InheritanceManager()

    email_weekly_updates = models.BooleanField(null=True, blank=True)
    email_discordance_updates = models.BooleanField(null=True, blank=True)
    columns = models.ForeignKey(CustomColumnsCollection, on_delete=SET_NULL, null=True, blank=True)
    default_sort_by_column = models.ForeignKey(CustomColumn, on_delete=SET_NULL, null=True, blank=True)
    variant_link_in_analysis_opens_new_tab = models.BooleanField(null=True)
    tool_tips = models.BooleanField(null=True, blank=True)
    node_debug_tab = models.BooleanField(null=True, blank=True)
    import_messages = models.BooleanField(null=True, blank=True)  # Get a message once import is done
    igv_port = models.IntegerField(null=True, blank=True)
    default_genome_build = models.ForeignKey(GenomeBuild, on_delete=CASCADE, null=True, blank=True)


class GlobalSettings(SettingsOverride):
    def save(self, *args, **kwargs):
        if not self.pk and GlobalSettings.objects.exists():
            raise ValidationError('There is can be only one GlobalSettings instance')
        return super(GlobalSettings, self).save(*args, **kwargs)

    def __str__(self):
        return "Global"


class OrganizationUserSettingsOverride(SettingsOverride):
    organization = models.OneToOneField(Organization, on_delete=CASCADE)

    def __str__(self):
        return str(self.organization)


class LabUserSettingsOverride(SettingsOverride):
    lab = models.OneToOneField(Lab, on_delete=CASCADE)

    def __str__(self):
        return str(self.lab)


class UserSettingsOverride(SettingsOverride):
    user = models.OneToOneField(User, on_delete=CASCADE)
    # records created by this user will also be associated to this lab by default
    default_lab = models.ForeignKey(Lab, null=True, blank=True, on_delete=SET_NULL)
    # Allows us to store OAuth sub against the user
    oauth_sub = models.TextField(null=True, blank=True)

    def auto_set_default_lab(self):
        user = self.user
        if self.default_lab:
            group_name = self.default_lab.group_name
            if not group_name:
                # default lab isn't associated with a group permission
                # so let it ride
                return

            # if user no longer belongs to the group associated to the lab
            # remove their default lab record
            if not user.groups.filter(name=group_name).exists():
                self.default_lab = None

        if not self.default_lab:
            try:
                self.default_lab = Lab.valid_labs_qs(self.user).first()
            except Lab.DoesNotExist:
                pass

    def __str__(self):
        return f"UserSettings for {self.user}"


@dataclass
class UserSettings:
    user: User
    email_weekly_updates: bool
    email_discordance_updates: bool
    columns: CustomColumnsCollection
    default_sort_by_column: CustomColumn
    variant_link_in_analysis_opens_new_tab: bool
    tool_tips: bool
    node_debug_tab: bool
    import_messages: bool
    igv_port: bool
    default_genome_build: GenomeBuild
    default_lab: Optional[Lab]
    oauth_sub: str
    _settings_overrides: List[SettingsOverride]

    @staticmethod
    def preferred_label_for(user: User) -> Optional[str]:
        if not user:
            return None
        if user.first_name and user.last_name:
            return user.first_name + ' ' + user.last_name
        return user.username

    @property
    def avatar_url(self) -> str:
        return avatar_url(self.user)

    @property
    def avatar_color(self) -> str:
        if self.user.username == 'admin_bot':
            return '#ffffff'
        colors = ['#ddddff', '#ffbbbb', '#aaddaa', '#dddddd', '#ddddbb']
        return colors[string_deterministic_hash(self.user.username) % len(colors)]

    @staticmethod
    def get_settings_overrides(user=None, lab=None, organization=None) -> List[SettingsOverride]:
        user_settings_override = None
        lab_settings_override = None

        if user:
            if not user.is_authenticated:
                raise ValueError("Cannot call UserSettings.get_for_user with unauth user!")
            user_settings_override = UserSettingsOverride.objects.get_or_create(user=user)[0]
            if lab is None:
                lab = user_settings_override.default_lab

        if lab:
            lab_settings_override = LabUserSettingsOverride.objects.get_or_create(lab=lab)[0]
            if organization is None:
                organization = lab.organization

        settings_overrides = [GlobalSettings.objects.get()]  # Later overrides earlier
        if organization:
            org_settings_override = OrganizationUserSettingsOverride.objects.get_or_create(organization=organization)[0]
            settings_overrides.append(org_settings_override)

        if lab_settings_override:
            settings_overrides.append(lab_settings_override)

        if user_settings_override:
            settings_overrides.append(user_settings_override)
        return settings_overrides

    @staticmethod
    def get_for_user(user: User) -> 'UserSettings':
        override_fields = [s.name for s in dataclasses.fields(UserSettings)]
        kwargs = {f: None for f in override_fields}  # Need to pass all params
        settings_overrides = UserSettings.get_settings_overrides(user=user)
        kwargs["_settings_overrides"] = settings_overrides
        for so in settings_overrides:
            for f in override_fields:
                val = getattr(so, f, None)
                if val is not None:
                    kwargs[f] = val
        return UserSettings(**kwargs)

    @lazy
    def initial_perm_read_and_write_groups(self) -> Tuple[Set[Group], Set[Group]]:
        groups = self.user.groups.all()
        settings_overrides = self._settings_overrides
        return self.get_initial_perm_read_and_write_groups(groups, settings_overrides)

    @staticmethod
    def get_initial_perm_read_and_write_groups(groups, settings_overrides) -> Tuple[Set[Group], Set[Group]]:
        group_read = defaultdict(lambda x: False)
        group_write = defaultdict(lambda x: False)
        qs = SettingsInitialGroupPermission.objects.filter(group__in=groups)
        for so in settings_overrides:
            for sigp in qs.filter(settings=so):
                if sigp.read is not None:
                    group_read[sigp.group] = sigp.read
                if sigp.write is not None:
                    group_write[sigp.group] = sigp.write

        read_groups = {g for g, read_perm in group_read.items() if read_perm}
        write_groups = {g for g, write_perm in group_write.items() if write_perm}
        return read_groups, write_groups

    def get_override_source_and_values_before_user(self) -> Tuple[Dict[str, str], Dict[str, str]]:
        override_fields = [s.name for s in dataclasses.fields(UserSettings)]
        parent_overrides = self._settings_overrides[:-1]  # Skip last, which is User override
        return self.get_override_source_and_values(override_fields, parent_overrides)

    @staticmethod
    def get_override_source_and_values(override_fields, parent_overrides) -> Tuple[Dict[str, str], Dict[str, str]]:
        override_source = {}
        override_values = {}
        for so in parent_overrides:
            source = str(so)
            for f in override_fields:
                val = getattr(so, f, None)
                if val is not None:
                    override_source[f] = source
                    override_values[f] = val
        return override_source, override_values

    def get_lab(self):
        """ If user has chosen a lab, use that.
            If they only have 1 use that """

        lab = self.default_lab
        if lab is None:
            valid_labs = list(Lab.valid_labs_qs(self.user))
            num_labs = len(valid_labs)
            if num_labs == 0:
                msg = "You do not belong to any labs."
                raise ValueError(msg)
            if num_labs == 1:
                lab = valid_labs[0]
            elif num_labs > 1:
                msg = f"You belong to {num_labs} labs. Default lab required."
                raise ValueError(msg)
        return lab

    def get_node_count_settings_collection(self) -> Optional['NodeCountSettingsCollection']:
        for so in reversed(self._settings_overrides):  # Just need last override
            try:
                return so.nodecountsettingscollection
            except NodeCountSettingsCollection.DoesNotExist:
                pass
        return None

    @staticmethod
    def get_genome_build_or_default(user: User, genome_build_name: Optional[str] = None) -> GenomeBuild:
        user_settings = UserSettings.get_for_user(user)
        if genome_build_name:
            genome_build = GenomeBuild.get_name_or_alias(genome_build_name)
        else:
            genome_build = user_settings.default_genome_build
        return genome_build

    @staticmethod
    def get_lab_and_error(user: User) -> Tuple[Optional[Lab], Optional[str]]:
        lab_error = None
        lab = None
        user_settings = UserSettings.get_for_user(user)
        try:
            lab = user_settings.get_lab()
        except ValueError as ve:
            lab_error = str(ve)
        return lab, lab_error

    @property
    def classification_issue_count(self) -> int:
        from classification.models import Classification
        from flags.models import FlagCollection
        return FlagCollection.filter_for_open_flags(
            Classification.filter_for_user(user=self.user)
        ).order_by('-created').exclude(withdrawn=True).count()

    def __str__(self):
        return f"UserSettings for {self.user}"


class SettingsInitialGroupPermission(models.Model):
    """ How to initially share objects w/Guardian permissions. """
    settings = models.ForeignKey(SettingsOverride, on_delete=CASCADE)
    group = models.ForeignKey(Group, on_delete=CASCADE)
    read = models.BooleanField(null=True, blank=True)
    write = models.BooleanField(null=True, blank=True)

    @staticmethod
    def create_global_settings(group: Group):
        """ Triggerd on Group post_save in app signal handler """
        global_settings = GlobalSettings.objects.get()
        SettingsInitialGroupPermission.objects.get_or_create(settings=global_settings,
                                                             group=group, read=True, write=False)


class NodeCountSettingsCollection(models.Model):
    # Need this object to distinguish empty counts vs not yet configured
    settings = models.OneToOneField(SettingsOverride, on_delete=CASCADE)

    def get_node_count_filters(self):
        qs = self.nodecountsettings_set.all().order_by("sort_order")
        return [nc.built_in_filter for nc in qs]


class AbstractNodeCountSettings(models.Model):
    built_in_filter = models.CharField(max_length=1, choices=BuiltInFilters.CHOICES, null=True)
    sort_order = models.IntegerField()

    class Meta:
        abstract = True

    @staticmethod
    def get_types_from_labels(node_count_labels):
        # Convert from labels to types
        labels = dict(BuiltInFilters.CHOICES)
        node_count_types = []
        for label in node_count_labels:
            type_info = {"label": labels[label]}
            if label == BuiltInFilters.TOTAL:
                type_info.update({"show_zero": True, "link": False})
            else:
                type_info.update({"show_zero": False, "link": True})

            node_count_types.append((label, type_info))
        return node_count_types

    @staticmethod
    def save_count_configs_from_array(record_set, node_counts_array):
        record_set.all().delete()  # Delete and recreate

        valid_choices = dict(BuiltInFilters.CHOICES)
        for i, nc in enumerate(node_counts_array):
            if nc in valid_choices:
                record_set.create(built_in_filter=nc,
                                  sort_order=i)


class NodeCountSettings(AbstractNodeCountSettings):
    collection = models.ForeignKey(NodeCountSettingsCollection, on_delete=CASCADE)


class UserContact(models.Model):
    """ Need to add phone - email is already built into User model """
    user = models.OneToOneField(User, on_delete=CASCADE)
    phone_number = models.TextField(null=True, blank=True)

    @staticmethod
    def get_for_user(user):
        user_contact, _ = UserContact.objects.get_or_create(user=user)
        return user_contact
