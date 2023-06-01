using Clapeyron

# Se incluye para que los alpha func de abajo se sumen a los de la Clapeyron
import Clapeyron: α_function

abstract type MathiasCopemanAlphaModel <: AlphaModel end

struct MathiasCopemanAlphaParam <: EoSParam
    c1::SingleParam{Float64}
    c2::SingleParam{Float64}
    c3::SingleParam{Float64}
end

@newmodelsimple MathiasCopemanAlpha MathiasCopemanAlphaModel MathiasCopemanAlphaParam
export MathiasCopemanAlpha

MathiasCopemanAlpha

function MathiasCopemanAlpha(components::Vector{String}; userlocations::Vector{String}=String[], verbose::Bool=false)
    params = getparams(components; userlocations=userlocations, verbose=verbose)
    c1 = params["c1"]
    c2 = params["c2"]
    c3 = params["c3"]
    packagedparams = MathiasCopemanAlphaParam(c1,c2,c3)
    model = MathiasCopemanAlpha(packagedparams, verbose=verbose)
    return model
end

MathiasCopemanAlpha

function α_function(model::CubicModel,V,T,z,alpha_model::MathiasCopemanAlphaModel)
    Tc = model.params.Tc.values
    _c1  = alpha_model.params.c1.values
    _c2  = alpha_model.params.c2.values
    _c3  = alpha_model.params.c3.values
    α = zeros(typeof(T*1.0),length(Tc))
    for i in Clapeyron.@comps
        c1 = _c1[i]
        c2 = _c2[i]
        c3 = _c3[i]
        Tr = T/Tc[i]

        if Tr < 1
            α[i] = (1 + c1 * (1 - (Tr)^0.5) + c2 * (1 - (Tr)^0.5)^2 + c3 * (1 - (Tr)^0.5)^3)^2
        else
            α[i] = (1 + c1 * (1 - (Tr)^0.5))^2
        end
    end
    return α
end

function α_function(model::CubicModel,V,T,z::Clapeyron.SingleComp,alpha_model::MathiasCopemanAlphaModel)
    Tc = model.params.Tc.values[1]
    c1  = alpha_model.params.c1.values[1]
    c2  = alpha_model.params.c2.values[1]
    c3  = alpha_model.params.c3.values[1]
    Tr = T/Tc

    if Tr < 1
        α = (1 + c1 * (1 - (Tr)^0.5) + c2 * (1 - (Tr)^0.5)^2 + c3 * (1 - (Tr)^0.5)^3)^2
    else
        α = (1 + c1 * (1 - (Tr)^0.5))^2
    end
end
