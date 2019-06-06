"""
    add_measure(model::InfiniteModel, v::Measure)
Add a measure to the `InfiniteModel` object in an analagous way to `JuMP.add_variable`.
"""
function add_measure(m::InfiniteModel, meas::Measure)
    m.next_meas_index += 1
    mref = MeasureRef(m, m.next_meas_index)
    m.measures[mref.index] = meas
    mref.model.meas_to_name[mref.index] = _make_meas_name(meas)
    return mref
end

# Set a default weight function
_w(t) = 1

# Define constructor functions
DiscretizeData() = DiscretizeData("measure", _w, [], [])
DiscretizeData(coeffs::Vector, pts::Array) = DiscretizeData("measure", _w, coeffs, pts)
DiscretizeData(name::String, coeffs::Vector, pts::Array) = DiscretizeData(name, _w, coeffs, pts)

"""
    measure(expr::Union{InfiniteExpr, MeasureRef}, data::AbstractMeasureData)
Implement a measure in an expression in a similar fashion to the `sum` method in JuMP.
"""
function measure(expr::Union{InfiniteExpr, MeasureRef}, data::AbstractMeasureData)
    meas = Measure(expr, data)
    model = _get_model_from_expr(expr)
    if model == nothing
        error("Expression contains no variables.")
    end
    return add_measure(model, meas)
end

# Parse the model pertaining to an expression
function _get_model_from_expr(expr::JuMP.AbstractJuMPScalar)
    if expr isa JuMP.AbstractVariableRef
        return expr.model
    elseif expr isa JuMP.GenericAffExpr
        aff_vars = [k for k in keys(expr.terms)]
        if length(aff_vars) > 0
            return aff_vars[1].model
        else
            return
        end
    elseif expr isa JuMP.GenericQuadExpr
        aff_vars = [k for k in keys(expr.aff.terms)]
        if length(aff_vars) > 0
            return aff_vars[1].model
        else
            var_pairs = [k for k in keys(expr.terms)]
            if length(var_pairs) > 0
                return var_pairs[1].a.model
            else
                return
            end
        end
    else
        return expr.m
    end
end

# Parse the string for displaying a measure
function _make_meas_name(meas::Measure)
    return string(meas.data.name, "(", JuMP.function_string(JuMP.REPLMode, meas.func), ")")
end

"""
    JuMP.name(mref::MeasureRef)
Extend the `JuMP.name` function to accomodate measure references.
"""
JuMP.name(mref::MeasureRef) = mref.model.meas_to_name[mref.index]

# TODO Add manipulation functions like variables have
