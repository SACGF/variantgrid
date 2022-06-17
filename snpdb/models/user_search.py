from typing import Any

from django.contrib.auth.models import User
from django.db.models import Q
from django.dispatch import receiver
from django.urls import reverse

from snpdb.search2 import SearchResponseRecordAbstract, search_signal, SearchInput, SearchResponse


class SearchResponseUser(SearchResponseRecordAbstract[User]):

    @classmethod
    def search_type(cls) -> str:
        return "Users"

    def get_absolute_url(self) -> str:
        return reverse('view_user', kwargs={"pk": self.record.pk})


@receiver(search_signal, sender=SearchInput)
def search_users(sender: Any, search_input: SearchInput, **kwargs) -> SearchResponse:
    response: SearchResponse[User] = SearchResponse(SearchResponseUser)
    if search_input.user.is_superuser and search_input.matches_has_alpha():
        response.extend(User.objects.filter(
            Q(username__icontains=search_input.search_string) or
            Q(first_name__icontains=search_input.search_string) or
            Q(last_name__icontains=search_input.search_string) or
            Q(email__icontains=search_input.search_string)
        ))
    return response
