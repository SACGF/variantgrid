from django.contrib.auth.models import User
from django.db.models import Q, Value, CharField
from django.db.models.functions import Lower, Concat
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
    user_combined_name_qs = User.objects.annotate(
        combined_name=Concat('first_name',  Value(' '), 'last_name', output_field=CharField())
    )
    for user in user_combined_name_qs.filter(
        search_input.q_words('username') |
        search_input.q_words('combined_name') |
        search_input.q_words('email')
    ).order_by('-is_active', Lower('username')):
        messages = []
        if not user.is_active:
            messages = ["User is inactive"]
        yield UserPreview(user), messages
