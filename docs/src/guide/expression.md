# Expressions
A guide for the defining and understanding the variable expressions
used in `InfiniteOpt`.

## Overview
Expressions in `InfiniteOpt` (also called functions) refer to mathematical
statements involving variables and numbers. Thus, these comprise the
mathematical expressions used that are used in measures, objectives, and
constraints. Programmatically, `InfiniteOpt` simply extends `JuMP` expression
types and methods principally pertaining to affine and quadratic mathematical
expressions. An natively supported abstraction for general nonlinear expressions
is currently under development since that of `JuMP` is not readily extendable.

## Datatype Hierarchy
Expressions employ variable reference types inherited from
[`JuMP.AbstractVariableRef`](@ref) to form expression objects. `InfiniteOpt`
uses a hierarchy of such types to organize the complexities associated with
modeling infinite dimensional programs. The figure below summarizes this
hierarchy of variable reference types where the abstract types are depicted in
green and the concrete types are shown blue.

![tree](../assets/variable_tree.png)

Following `JuMP`, expression objects are parameterized by the variable reference
type that is present in the expression. In `InfiniteOpt` expressions
automatically, select the most specific variable reference type possible in
accordance with the above figure. For instance, an expression that only contains
hold variables will be classified as a [`HoldVariableRef`](@ref) expression object,
whereas an expression containing hold variables and a measure would be classified
as a [`MeasureFiniteVariableRef`](@ref) expression object. This hierarchical
classification becomes convenient to guide infinite program reformulation
schemes in how to treat different expressions. The default transcription
methodology employed by `InfiniteOpt.TranscriptionOpt` uses these classifications
to efficiently differentiate between finite and infinite expressions.

## Affine Expressions
An affine expression pertains to a mathematical function of the form:
```math
f_a(x) = a_1x_1 + ... + a_nx_n + b
```
where ``x \in \mathbb{R}^n`` denote variables, ``a \in \mathbb{R}^n`` denote
coefficients, and ``b \in \mathbb{R}`` denotes a constant value. Such
expressions, are prevalent in any problem than involves linear constraints
and/or objectives.

In `InfiniteOpt`, affine expressions can be defined directly
using `Julia`'s arithmetic operators (i.e., `+`, `-`, `*`, etc.) or using
[`@expression`](@ref).  

To be continued...

## Quadratic Expressions


## Nonlinear Expressions
