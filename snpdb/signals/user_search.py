from functools import cached_property
from typing import Any, Type, Optional, List
from django.contrib.auth.models import User
from django.db.models import Q
from django.dispatch import receiver
from library.preview_request import PreviewData
from snpdb.models import AvatarDetails
from snpdb.search2 import SearchResponseRecordAbstract, search_signal, SearchInput, SearchResponse


class SearchResponseUser(SearchResponseRecordAbstract[User]):

    @classmethod
    def result_class(cls) -> Type:
        return User

    @classmethod
    def category(cls) -> str:
        return "User"

    @cached_property
    def preview(self) -> PreviewData:
        return AvatarDetails.avatar_for(self.record).preview()

    @property
    def messages(self) -> Optional[List[str]]:
        if not self.record.is_active:
            return ["User is inactive"]


@receiver(search_signal, sender=SearchInput)
def user_search(sender: Any, search_input: SearchInput, **kwargs) -> SearchResponse:
    response: SearchResponse[User] = SearchResponse(SearchResponseUser)
    if search_input.user.is_superuser and search_input.matches_has_alpha():
        response.extend(User.objects.filter(
            Q(username__icontains=search_input.search_string) or
            Q(first_name__icontains=search_input.search_string) or
            Q(last_name__icontains=search_input.search_string) or
            Q(email__icontains=search_input.search_string)
        ))
    return response
