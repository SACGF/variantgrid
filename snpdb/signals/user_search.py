from django.contrib.auth.models import User
from django.db.models import Q
from snpdb.models import UserPreview
from snpdb.search import search_receiver, SearchInputInstance, \
    HAS_ALPHA_PATTERN, SearchExample


#FIXME make a new class for search result for user that can be previewable

@search_receiver(
    search_type=UserPreview,
    pattern=HAS_ALPHA_PATTERN,
    admin_only=True,
    example=SearchExample(
        note="Search on username, name or email",
        examples=["jane@institute.org.au", "Lisa"]
    )
)
def user_search(search_input: SearchInputInstance):
    # TODO search for search-words across all fields (so you can search for first name last name for example)
    for user in User.objects.filter(
        Q(username__icontains=search_input.search_string) or
        Q(first_name__icontains=search_input.search_string) or
        Q(last_name__icontains=search_input.search_string) or
        Q(email__icontains=search_input.search_string)
    ):
        messages = []
        if not user.is_active:
            messages = ["User is inactive"]
        yield UserPreview(user), messages
