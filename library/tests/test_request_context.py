from django.http import StreamingHttpResponse
from django.test import RequestFactory, SimpleTestCase

from library.request_context import (
    RequestContextMiddleware,
    get_current_request,
    get_request_variable,
    set_thread_variable,
)


class RequestContextMiddlewareTest(SimpleTestCase):

    def tearDown(self):
        set_thread_variable("request", None)

    def test_streamed_content_sees_the_request_until_the_response_closes(self):
        """ Streaming exports read the current request while the content is generated, after the view returned """
        def get_response(_request):
            return StreamingHttpResponse(str(get_current_request() is not None) for _ in range(1))

        request = RequestFactory().get("/")
        response = RequestContextMiddleware(get_response)(request)
        self.assertEqual(b"".join(response.streaming_content), b"True")
        self.assertIs(get_current_request(), request)

        response.close()
        self.assertIsNone(get_current_request())

    def test_request_only_cache_refuses_without_a_request(self):
        """ QuerySetRequestCache relies on this so a celery worker never caches across tasks """
        with self.assertRaises(RuntimeError):
            get_request_variable("key", use_threadlocal_if_no_request=False)
