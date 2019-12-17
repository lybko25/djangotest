# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from rest_framework.views import APIView
from rest_framework.response import Response
import json
from django.core import serializers
import  numpy as np
import sympy as sym
from django.shortcuts import render
from sympy.parsing.sympy_parser import parse_expr
# Create your views here.

from .models import Equation
eps = 10**-4;
class EquationView(APIView):
    def get(self, request):
        equations = serializers.serialize("json",Equation.objects.all(),fields=('id','equation', 'result'))
        return Response({"equations": equations})
    def post(self, request):
        
        x1, x2,x3 = sym.symbols('x1 x2 x3')
        eps = request.data.get('eps')
        method = request.data.get('method')
        symbols1 = [x1,x2]
        f1 = sym.sympify(request.data.get('function'))
        f2= [];
        g1 = sym.sympify(request.data.get('gFunction'))
        g2=  []
        h2 = []
        x0 = sym.sympify(request.data.get('x0'))
        xm1 = sym.sympify(request.data.get('xm1'))
        x01 = []
        xm11 =[]
        for x in range(len(f1)):
            f2.append(parse_expr(f1[x]))
            g2.append(parse_expr(g1[x]))
            h2.append(f2[x] + g2[x])
            x01.append(x0[x])
            xm11.append(xm1[x])

        f1=sym.Array(f2)
        g1=sym.Array(g2)
        h1=sym.Array(h2)
        # return Response({"dsadas":  "json", "result": str(f1),"result1": str(g1),"result2": str(h1)})

        x0=sym.Matrix(x01)
        xm1=sym.Matrix(xm11)
        result = 0;
        if method == 'oneStep':
            result = twoStepMethod2(f1,g1,h1, symbols1,x0,xm1)
        elif method == 'twoStep':
            result = twoStepMethod1(f1,g1,h1, symbols1,x0,xm1)
        elif method == 'chords':
            result = twoStepMethod3(f1,g1,h1, symbols1,x0,xm1)
        elif method == 'kurchatova':
            result = twoStepMethod4(f1,g1,h1, symbols1,x0,xm1)
        return Response({"dsadas":  "json", "result": str(result)})


def Jacobian(v_str, f_list, x):
    array = {}
    for a in range(len(x)):
        array["x"+str(a+1)] = x[a];
    
    J = sym.zeros(len(f_list),len(v_str))
    for i, fi in enumerate(f_list):
        for j, s in enumerate(v_str):
            J[i,j] = sym.diff(fi, s).evalf(subs=array)
    return J

def dividedDifferences(x,y, gFuncs):
    n = len(gFuncs);
    gMatrix = sym.zeros(n);
    for i in range(n):
        for j in range(n):
            subs1 = {};            
            subs2 = {};

            for a in range(n):
                if j <= a:
                    subs1["x"+str(a+1)] = x[a];
                elif j > a:
                    subs1["x"+str(a+1)] = y[a];
                if j<a:
                    subs2["x"+str(a+1)] = x[a];
                elif j>=a:
                    subs2["x"+str(a+1)] = y[a];

            gMatrix[i,j]= (gFuncs[i].evalf(subs = subs1)-gFuncs[i].evalf(subs=subs2))/(x[j]-y[j]);
    
    return gMatrix;

def checkResult(hFunc, x0):
    subs = {}
    for a in range(len(x0)):
        subs["x"+str(a+1)] = x0[a];
    
    for i in range(len(hFunc)):
        if sym.Abs(hFunc[i].evalf(subs = subs)) > eps:
            return True
    return False
        
def calculateFunc(funcs, subs):
    resArray = sym.zeros(len(funcs), 1);
    
    for i in range(len(funcs)):
        resArray[i] = funcs[i].evalf(subs = subs);
    
    return resArray;
        
    
#two-step
def twoStepMethod1(f, g, h, s, x0, y0):
    x = x0;
    y = y0;
    iterations = 0;
    while(checkResult(h, x)):
        iterations = iterations + 1
        subX = {}        
        subY = {}

        for a in range(len(x)):
            subX["x"+str(a+1)] = x[a];
        An = (Jacobian(s, f, (x + y)/2) + dividedDifferences(x,y,g));
        x = x - An.inv()*(calculateFunc(f,subX) + calculateFunc(g, subX));
        for a in range(len(x)):
            subY["x"+str(a+1)] = x[a];
        y = x-An.inv()*(calculateFunc(f,subY) + calculateFunc(g, subY))
    return x;

#oneStep
def twoStepMethod2(f, g, h, s, x0, y0):
    x = x0;
    y = y0;
    iterations = 0;
    while(checkResult(h, x)):
        iterations = iterations + 1
        subX = {}        

        for a in range(len(x)):
            subX["x"+str(a+1)] = x[a];
        An = Jacobian(s, f, x);
        x = x - An.inv()*(calculateFunc(f,subX) + calculateFunc(g, subX)); 
    print(x,' iterations = ', iterations)
    print()


#chords 
def twoStepMethod3(f, g, h, s, x0, y0):
    x = x0;
    xn1 = y0;
    iterations = 0;
    while(checkResult(h, x)):
        iterations = iterations + 1
        subX = {}        
        subY = {}

        for a in range(len(x)):
            subX["x"+str(a+1)] = x[a];
        An = (Jacobian(s, f, x) + dividedDifferences(x,xn1,g));
        xn1 = x
        x = x - An.inv()*(calculateFunc(f,subX) + calculateFunc(g, subX));
        
        for a in range(len(x)):
            subY["x"+str(a+1)] = x[a];
    print(x,' iterations = ', iterations)
    print()
    
#kurchatova
def twoStepMethod4(f, g, h, s, x0, y0):
    x = x0;
    xn1 = y0;
    iterations = 0;
    while(checkResult(h, x)):
        iterations = iterations + 1
        subX = {}        
        subY = {}

        for a in range(len(x)):
            subX["x"+str(a+1)] = x[a];
        An = (Jacobian(s, f, x) + dividedDifferences(2*x-xn1,xn1,g));
        xn1 = x
        x = x - An.inv()*(calculateFunc(f,subX) + calculateFunc(g, subX));
        
        for a in range(len(x)):
            subY["x"+str(a+1)] = x[a];
    print(x,' iterations = ', iterations)
    print()
