from typing import Any
from django.contrib.auth.models import User
from django.db.models import Q
from django.dispatch import receiver
from snpdb.models import AvatarDetails
from snpdb.search2 import search_signal, SearchInput, SearchResponse


@receiver(search_signal, sender=SearchInput)
def user_search(sender: Any, search_input: SearchInput, **kwargs) -> SearchResponse:
    if search_input.user.is_superuser and search_input.matches_has_alpha():
        response = SearchResponse("User")
        for user in User.objects.filter(
            Q(username__icontains=search_input.search_string) or
            Q(first_name__icontains=search_input.search_string) or
            Q(last_name__icontains=search_input.search_string) or
            Q(email__icontains=search_input.search_string)
        ):
            messages = []
            if not user.is_active:
                messages = ["User is inactive"]
            response.add(AvatarDetails.avatar_for(user).preview(), messages=messages)

        return response
