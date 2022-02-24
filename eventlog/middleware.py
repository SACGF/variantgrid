class PageViewsMiddleware:

    def process_view(self, request, view_func, view_args, view_kwargs):
        if name := view_func.__name__:
            print(f"{request.user} just views {name} with args {view_args}")

        pass