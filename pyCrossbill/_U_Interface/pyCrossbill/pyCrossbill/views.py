from django.shortcuts import render
from django.http import HttpResponse

# create views here:

def first(request):
    return render(request, "first.html")