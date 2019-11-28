from django.urls import path
from .views import EquationView

app_name = "equations"
# app_name will help us do a reverse look-up latter.
urlpatterns = [
    path('equations/', EquationView.as_view()),
]