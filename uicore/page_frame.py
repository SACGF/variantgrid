"""
The page frame (#2007): everything around a page's content that depends on who is asking - the top bar and side bar,
username, avatar title, inbox count, site messages, Django messages and the Rollbar person. A page body carries none
of it; global.js fetches it from uicore/views/page_frame_view.py:page_frame after the page loads, so the same body can
serve every user (#695).

The menus depend only on the url name and the MenuRole, never on the individual user, so menu_html is memoised per
process (the registry and settings change only on restart). The rest is per user and built per request.

Entry points: get_page_frame(request, url_name), menu_html(url_name, role).
"""
from dataclasses import dataclass, field
from enum import Enum
from functools import lru_cache
from typing import Optional

from django.conf import settings
from django.contrib.messages import get_messages
from django.http import HttpRequest
from django.template.loader import render_to_string
from django.urls import reverse

from snpdb.models.models import SiteMessage
from snpdb.models.models_user_settings import AvatarDetails
from snpdb.user_settings_manager import UserSettingsManager
from uicore.menus import MenuItem, current_menu, get_menus
from user_messages.models import inbox_count_for
from variantgrid.perm_path import get_visible_url_names


class MenuRole(Enum):
    """ Everything the menus depend on about the user: admin-only items need a superuser """
    ANONYMOUS = 'anonymous'
    USER = 'user'
    SUPERUSER = 'superuser'

    @staticmethod
    def for_user(user) -> 'MenuRole':
        if user.is_superuser:
            return MenuRole.SUPERUSER
        if user.is_authenticated:
            return MenuRole.USER
        return MenuRole.ANONYMOUS


@dataclass(frozen=True)
class MenuHtml:
    main: str = ''  # top bar
    sub: str = ''  # side bar for the page


@dataclass
class PageFrame:
    menu_main_html: str
    menu_sub_html: str
    user_html: str  # inbox link and username / avatar title, '' for anonymous
    username: str  # '' for anonymous
    site_messages_html: str
    messages_html: str  # Django messages, consumed by this request
    rollbar_person: dict = field(default_factory=dict)  # {} for anonymous


def _item_context(item: MenuItem, active: bool) -> dict:
    return {
        'url': item.href or reverse(item.url_name),
        'css_class': item.css_class,
        'title': item.display_title,
        'icon': item.icon,
        'admin_only': item.admin_only,
        'active': active,
        'type': 'side',
        'id': f'submenu-{item.url_name}',
        'external': item.external,
        'method': item.method,
    }


def _menu_bar_main(active_menu) -> str:
    top_items = [{
        'url': reverse(menu.url_name),
        'title': menu.title,
        'active': menu == active_menu,
        'type': 'top',
        'id': f'menu-top-{menu.key}',
        'method': 'get',
    } for menu in get_menus() if menu.in_top_bar]
    return render_to_string("uicore/menus/menu_bar_main.html", {
        'top_items': top_items,
        'help_url': settings.HELP_URL,
    })


def _menu_bar_sub(menu, url_name: str, role: MenuRole) -> str:
    if not menu:
        return ''
    items = [_item_context(item, item.owns(url_name)) for item in menu.visible_items()
             if role == MenuRole.SUPERUSER or not item.admin_only]
    return render_to_string("uicore/menus/menu_bar_sub.html", {
        'items': items,
        'footer_template': menu.footer_template,
    })


@lru_cache(maxsize=2048)
def menu_html(url_name: Optional[str], role: MenuRole) -> MenuHtml:
    """ Anonymous users only reach base_external pages, so they get no menus """
    if role == MenuRole.ANONYMOUS:
        return MenuHtml()
    if url_name not in get_visible_url_names():
        url_name = None  # keeps the memo bounded whatever the client sends
    menu = current_menu(url_name)
    return MenuHtml(main=_menu_bar_main(menu), sub=_menu_bar_sub(menu, url_name, role))


def _messages_html(request: HttpRequest) -> str:
    if messages := list(get_messages(request)):
        return render_to_string("uicore/messages/messages.html", {'messages': messages})
    return ''


def _user_html(user) -> str:
    avatar_details = AvatarDetails.avatar_for(user)
    title_icon_html = ''
    # titles first: most users hold none, and their UserSettings costs several queries
    if avatar_details.titles and avatar_details.shows_titles_for(UserSettingsManager.get_user_settings(user)):
        title_icon_html = avatar_details.title_icon_html
    return render_to_string("uicore/page/navbar_user.html", {
        'user': user,
        'inbox_enabled': settings.INBOX_ENABLED,
        'mail_count': inbox_count_for(user) if settings.INBOX_ENABLED else 0,
        'title_icon_html': title_icon_html,
    })


def get_page_frame(request: HttpRequest, url_name: Optional[str]) -> PageFrame:
    user = request.user
    menus = menu_html(url_name, MenuRole.for_user(user))
    frame = PageFrame(menu_main_html=menus.main, menu_sub_html=menus.sub, user_html='', username='',
                        site_messages_html='', messages_html=_messages_html(request))
    if user.is_authenticated:
        frame.user_html = _user_html(user)
        frame.username = user.username
        frame.site_messages_html = render_to_string("uicore/site_messages/site_messages.html",
                                                     {'site_messages': SiteMessage.get_site_messages()})
        frame.rollbar_person = {'id': user.pk, 'username': user.username, 'email': user.email}
    return frame
