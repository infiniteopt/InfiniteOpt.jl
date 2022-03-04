################################################################################
#                       BASIC DERIVATIVE EVALUATION TYPES
################################################################################
"""
    AbstractDerivativeMethod

An abstract type for storing derivative evaluation data that is pertinent to its 
reformation/transcription. 
"""
abstract type AbstractDerivativeMethod end 

"""
    GenerativeDerivativeMethod <: AbstractDerivativeMethod

An abstract type for derivative evaluation method types that will require support 
generation when employed (e.g., internal node points associated with orthogonal 
collocation). Such methods can be used with derivatives that depend on independent 
infinite parameters, but cannot be used for ones that depend on dependent parameters.
"""
abstract type GenerativeDerivativeMethod <: AbstractDerivativeMethod end 

"""
    NonGenerativeDerivativeMethod <: AbstractDerivativeMethod

An abstract type for derivative evaluation method types that do not require the 
definition of additional support points. Such methods are amendable to any 
derivative in InfiniteOpt including those with dependent infinite parameter 
dependencies.
"""
abstract type NonGenerativeDerivativeMethod <: AbstractDerivativeMethod end

"""
    FDTechnique

An abstract data type for labels of specific techniques applied in the finite 
difference method in derivative evaluation.
"""
abstract type FDTechnique end

"""
    Forward <: FDTechnique

A technique label for finite difference method that implements a forward 
difference approximation.
"""
struct Forward <: FDTechnique end

"""
    Central <: FDTechnique

A technique label for finite difference method that implements a central 
difference approximation.
"""
struct Central <: FDTechnique end

"""
    Backward <: FDTechnique

A technique label for finite difference method that implements a backward 
difference approximation.
"""
struct Backward <: FDTechnique end

"""
    FiniteDifference{T <: FDTechnique} <: NonGenerativeDerivativeMethod

A `DataType` for information about finite difference method applied to 
a derivative evaluation. Note that the constructor is of the form:
```julia 
    FiniteDifference([technique::FDTechnique = Backward()],
                     [add_boundary_constr::Bool = true])
```
where `technique` is the indicated finite difference method to be applied and 
`add_boundary_constr` indicates if the finite difference equation corresponding to 
a boundary support should be included. Thus, for backward difference since
corresponds to the terminal point and for forward difference this corresponds to 
the initial point. We recommend using `add_boundary_constr = false` when an final 
condition is given with a backward method or when an initial condition is given 
with a forward method. Note that this argument is ignored for central finite 
difference which cannot include any boundary points.

**Fields** 
- `technique::T`: Mathematical technqiue behind finite difference
- `add_boundary_constraint::Bool`: Indicate if the boundary constraint should be 
  included in the transcription (e.g., the terminal boundary backward equation for 
  backward difference)
"""
struct FiniteDifference{T <: FDTechnique} <: NonGenerativeDerivativeMethod 
    technique::T
    add_boundary_constraint::Bool
    # set the constructor 
    function FiniteDifference(technique::T = Backward(), 
        add_boundary_constr::Bool = true) where {T <: FDTechnique}
        return new{T}(technique, add_boundary_constr)
    end
end

################################################################################
#                        DERIVATIVE METHOD FUNCTIONS
################################################################################
# Determine if any derivatives have derivative constraints
function has_derivative_constraints(pref::IndependentParameterRef)::Bool
    return _data_object(pref).has_deriv_constrs
end

# Make update function for whether it has derivative supports 
function _set_has_derivative_constraints(pref::IndependentParameterRef, 
                                         status::Bool)::Nothing 
    _data_object(pref).has_deriv_constrs = status
    return
end

"""
    derivative_method(pref::IndependentParameterRef)::AbstractDerivativeMethod

Returns the numerical derivative evaluation method employed with `pref` when it 
is used as an operator parameter in a derivative.

**Example**
```julia-repl
julia> derivative_method(pref) 
FiniteDifference(Backward, true)
```
"""
function derivative_method(pref::IndependentParameterRef)::AbstractDerivativeMethod
    return _core_variable_object(pref).derivative_method
end

# Make method to reset derivative constraints (supports are handled separately)
function _reset_derivative_constraints(pref::Union{IndependentParameterRef, 
                                                   DependentParameterRef})::Nothing
    if has_derivative_constraints(pref)
        @warn("Support/method changes will invalidate existing derivative evaluation " *
              "constraints that have been added to the InfiniteModel. Thus, " *
              "these are being deleted.")
        for idx in _derivative_dependencies(pref)
            delete_derivative_constraints(DerivativeRef(JuMP.owner_model(pref), idx))
        end
        _set_has_derivative_constraints(pref, false)
    end
    return
end

"""
    set_derivative_method(pref::IndependentParameterRef, 
                          method::AbstractDerivativeMethod)::Nothing

Specfies the desired derivative evaluation method `method` for derivatives that are 
taken with respect to `pref`. Any internal supports exclusively associated with 
the previous method will be deleted. Also, if any derivatives were evaluated 
manually, the associated derivative evaluation constraints will be deleted. Errors 
if new derivative method generates supports that are incompatible with existing 
measures.

**Example**
```julia-repl
julia> set_derivative_method(d, OrthogonalCollocation(2))

```
"""
function set_derivative_method(pref::IndependentParameterRef, 
    method::NonGenerativeDerivativeMethod
    )::Nothing
    old_param = _core_variable_object(pref)
    domain = _parameter_domain(pref)
    supps = _parameter_supports(pref)
    sig_figs = significant_digits(pref)
    if isempty(_generative_measures(pref))
        _reset_generative_supports(pref)
        new_param = IndependentParameter(domain, supps, sig_figs, method, 
                                         NoGenerativeSupports())
    else
        info = generative_support_info(pref)
        new_param = IndependentParameter(domain, supps, sig_figs, method, info)
    end
    _reset_derivative_constraints(pref)
    _set_core_variable_object(pref, new_param)
    if is_used(pref)
        set_optimizer_model_ready(JuMP.owner_model(pref), false)
    end
    return
end

# GenerativeDerivativeMethod
function set_derivative_method(pref::IndependentParameterRef, 
    method::GenerativeDerivativeMethod
    )::Nothing
    new_info = generative_support_info(method)
    old_info = generative_support_info(pref)
    if !isempty(_generative_measures(pref)) && new_info != old_info 
        error("Generative derivative method conflicts with existing generative " *
              "measures.")
    end
    old_param = _core_variable_object(pref)
    domain = _parameter_domain(pref)
    supps = _parameter_supports(pref)
    sig_figs = significant_digits(pref)
    new_param = IndependentParameter(domain, supps, sig_figs, method, new_info)
    _reset_derivative_constraints(pref)
    _reset_generative_supports(pref)
    _set_core_variable_object(pref, new_param)
    if is_used(pref)
        set_optimizer_model_ready(JuMP.owner_model(pref), false)
    end
    return
end