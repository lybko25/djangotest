# -*- coding: utf-8 -*-
from __future__ import unicode_literals

from django.db import models

# Create your models here.
class Equation(models.Model):
  id = models.IntegerField(primary_key=True)
  equation = models.CharField(max_length=1000)
  result = models.CharField(max_length=250)