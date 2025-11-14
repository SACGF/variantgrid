from django.shortcuts import render


def reanalyis(request):
    context = {}
    return render(request, 'analysis/reanalysis/reanalysis.html', context)


def view_reanalyis(request):
    context = {}
    return render(request, 'analysis/reanalysis/view_reanalysis.html', context)


