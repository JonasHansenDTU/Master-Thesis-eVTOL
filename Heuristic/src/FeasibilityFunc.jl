mutable struct planeSolution
    flightLegs::Int32

    route::Array{Int32,1}

    turnaroundTime::Array{Int32,1}
end

mutable struct allPlaneSolution
    planes::Vector{planeSolution}
end

function BatteryNeeded(TravelLength::Float32, battery_per_km::Float32)
    return TravelLength * battery_per_km
end

function BatteryCharged(turnaroundTime::Float16, ec::Float32)
    return turnaroundTime*ec
end

function FeasibleBattery(evtols::allPlaneSolution, bmax::Float32, bmid::Float32, bmin::Float32, dist::Dict{Tuple{Int,Int},Float64}, ec::Float32, battery_per_km::Float32, b_penalty)

    battery_overrule = 0

    for evtol in evtols.planes
        BatteryLevel = zeros(Float32, evtol.flightLegs + 1)
        BatteryLevel[1] = bmid

        for i in 1:evtol.flightLegs
            from = evtol.route[i]
            to = evtol.route[i + 1]
            TravelLength = Float32(dist[(from, to)])

            BatteryLevel[i + 1] =
                min(BatteryLevel[i] + BatteryCharged(Float16(evtol.turnaroundTime[i]), ec), bmax) -
                BatteryNeeded(TravelLength, battery_per_km)

            if BatteryLevel[i+1] > 80
                for i in bmid+1:BatteryLevel[i+1]
                    battery_overrule += b_penalty
                end
            end

            if BatteryLevel[i + 1] < bmin
                return false
            end
        end
    end 

    return true, battery_overrule
end

function FeasibleCompletionTime(evtols::allPlaneSolution, rt::Matrix{Int}, ET::Int)
    for evtol in evtols.planes
        travelTime = 0

        for i in 1:evtol.flightLegs
            from = evtol.route[i]
            to = evtol.route[i+1]

            travelTime += evtol.turnaroundTime[i] + rt[from, to]

            if travelTime > ET
                return false
            end
        end
    end

    return true
end

function FeasibleVertiportCapacity(evtols::allPlaneSolution, rt::Matrix{Int}, T::Int, V::Int, cap_v::Dict{Int64, Int64}, ET::Int)
    parkingTimes = zeros(Int, V, T)

    for evtol in evtols.planes
        travelTime = 0

        for i in 1:evtol.flightLegs
            from = evtol.route[i]
            to = evtol.route[i+1]

            startTime = travelTime
            endTime = travelTime + evtol.turnaroundTime[i]

            if endTime > ET 
                endTime = ET
            end

            for t in startTime:endTime
                if t >= 1
                    parkingTimes[from, t] += 1
                    if parkingTimes[from, t] > cap_v[from]
                        return false
                    end
                end
            end

            travelTime += evtol.turnaroundTime[i] + rt[from, to]
        end
    end

    return true
end

# function FeasibleCorridor(evtols::allPlaneSolution, rt::Matrix{Int}, T::Int, V::Int, cap_flt::Int, ET::Int)
#     destinationTimes = zeros(Int, V, V, T)

#     for evtol in evtols.planes
#         travelTime = 0

#         for i in 1:evtol.flightLegs
#             from = evtol.route[i]
#             to = evtol.route[i+1]

#             startTime = travelTime + evtol.turnaroundTime[i]
#             endTime = travelTime + evtol.turnaroundTime[i] + rt[from, to]

#             if endTime > ET 
#                 endTime = ET
#             end

#             for t in startTime:endTime
#                 if t >= 1
#                     destinationTimes[from, to, t] += 1
#                     if destinationTimes[from, to, t] > cap_flt 
#                         return false
#                     end
#                 end
#             end

#             travelTime += evtol.turnaroundTime[i] + rt[from, to]
#         end
#     end

#     return true
# end

function FeasibilityCheck(bmax::Float32, bmid::Float32, bmin::Float32, 
    dist::Dict{Tuple{Int,Int},Float64},
    ec::Float32, battery_per_km::Float32,
    evtols::allPlaneSolution,
    rt::Matrix{Int}, ET::Int, T::Int, V::Int,
    # cap_flt::Int, 
    cap_v::Dict{},
    b_penalty)

    P = zeros(Int32, 4)

    if FeasibleBattery(evtols, bmax, bmid, bmin, dist, ec, battery_per_km, b_penalty) == false
        P[1] = 1_000_000
    end

    if FeasibleBattery(evtols, bmax, bmid, bmin, dist, ec, battery_per_km, b_penalty) == true
        battery_overrule = FeasibleBattery(evtols, bmax, bmid, bmin, dist, ec, battery_per_km, b_penalty)
        P[1] = battery_overrule
    end

    if FeasibleCompletionTime(evtols, rt, ET) == false
        P[2] = 1_000_000
    end

    if FeasibleVertiportCapacity(evtols, rt, T, V, cap_v, ET) == false
        P[3] = 1_000_000
    end

    # if FeasibleCorridor(evtols, rt, T, V, cap_flt, ET) == false
    #     P[4] = 1
    # end

    return P
end