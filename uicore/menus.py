"""
Menus as data (#2007): Menu and MenuItem describe the top bar and every sub-menu, and a page's sub-menu and
highlighted items follow from its url name - templates never choose a menu. The declaration itself is the
deployment's `settings.MENUS` (an import path; the default is variantgrid/menus.py:MENUS), so a deployment can
extend or replace the registry without touching this app.

A url name belongs to a menu as an item (listed in the sub-menu), as one of an item's `pages` (a detail page that
highlights that item, e.g. view_liftover_run under Liftover) or as one of the menu's own `pages` (in the menu but
under no item, e.g. view_allele under Variants). Per-deployment differences stay in URLS_NAME_REGISTER
(variantgrid/perm_path.py:get_visible_url_names); `condition` is only for the few items that depend on a setting.

Entry points: get_menus(), current_menu(url_name); uicore/page_frame.py:menu_html renders them.
"""
from collections.abc import Callable
from dataclasses import dataclass
from typing import Optional

from django.conf import settings
from django.utils.module_loading import import_string

from variantgrid.perm_path import get_visible_url_names


def is_url_visible(url_name: str) -> bool:
    """ Registered on this deployment (settings.URLS_NAME_REGISTER) """
    return get_visible_url_names()[url_name]


@dataclass(frozen=True)
class MenuItem:
    url_name: str
    title: Optional[str] = None  # defaults to the url name in title case
    admin_only: bool = False
    icon: Optional[str] = None
    href: Optional[str] = None  # link to this rather than reversing url_name
    external: bool = False
    method: str = 'get'  # 'post' renders a hidden form, for views that change state
    css_class: str = ''
    pages: tuple[str, ...] = ()  # url names of detail pages that highlight this item
    condition: Optional[Callable[[], bool]] = None  # the item isn't in this menu at all when False

    @property
    def is_available(self) -> bool:
        return self.condition is None or self.condition()

    @property
    def is_visible(self) -> bool:
        return self.is_available and (bool(self.href) or is_url_visible(self.url_name))

    @property
    def display_title(self) -> str:
        return self.title or self.url_name.replace('_', ' ').title()

    def owns(self, url_name: str) -> bool:
        return url_name == self.url_name or url_name in self.pages


@dataclass(frozen=True)
class Menu:
    key: str
    title: str
    url_name: Optional[str]  # the top bar link; None keeps the menu out of the top bar (Settings)
    items: tuple[MenuItem, ...]
    pages: tuple[str, ...] = ()  # url names shown with this sub-menu without highlighting an item
    top_bar_condition: Optional[Callable[[], bool]] = None
    footer_template: Optional[str] = None

    @property
    def in_top_bar(self) -> bool:
        if not self.url_name or (self.top_bar_condition and not self.top_bar_condition()):
            return False
        return is_url_visible(self.url_name)

    def visible_items(self) -> list[MenuItem]:
        return [item for item in self.items if item.is_visible]

    def owns(self, url_name: str) -> bool:
        """ Ignores URLS_NAME_REGISTER: a superuser reaching a hidden page still gets its menu """
        return url_name in self.pages or any(item.owns(url_name) for item in self.items if item.is_available)


def get_menus() -> tuple[Menu, ...]:
    return import_string(settings.MENUS)


def current_menu(url_name: Optional[str]) -> Optional[Menu]:
    """ The menu that owns url_name, preferring one that is reachable (in the top bar, or has no top bar entry) -
        only Liftover and Seq / Software Versions are in two menus, and move to Settings when their own is off """
    if not url_name:
        return None
    owners = [menu for menu in get_menus() if menu.owns(url_name)]
    for menu in owners:
        if menu.url_name is None or menu.in_top_bar:
            return menu
    return next(iter(owners), None)
