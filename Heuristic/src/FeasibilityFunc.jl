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
    battery_levels = Vector{Vector{Float32}}()
    battery_overrule = 0

    feasible = true

    for evtol in evtols.planes
        BatteryLevel = zeros(Float32, evtol.flightLegs + 1)
        BatteryLevel[1] = bmid

        if evtol.flightLegs >0
            from2 = evtol.route[1]
            to2 = evtol.route[2]
            TravelLength2 = Float32(dist[(from2, to2)])

            # First leg: eVTOL starts at bmid. It charges during the initial
            # turnaround (ground time at base before first departure) ONLY if the
            # current battery is insufficient to fly the leg and still land at or
            # above bmin. The amount charged is the shortfall relative to the
            # usable battery above bmin, capped by available time and headroom.
            need2     = BatteryNeeded(TravelLength2, battery_per_km)
            shortfall = need2 - (BatteryLevel[1] - bmin)
            time_poss = BatteryCharged(Float16(evtol.turnaroundTime[1]), ec)
            headroom  = bmax - BatteryLevel[1]
            charge2   = min(max(0.0f0, shortfall), time_poss, headroom)

            BatteryLevel[2] = BatteryLevel[1] + charge2 - need2

            # Penalty if the peak level (after charging, before flying) exceeds bmid
            if BatteryLevel[2] + need2 > bmid
                for j in bmid+1:(BatteryLevel[2] + need2)
                    battery_overrule += b_penalty
                end
            end

            # Feasibility: post-flight level must stay at or above bmin
            if BatteryLevel[2] < bmin
                feasible = false
            end
        end

        for i in 2:evtol.flightLegs
            from = evtol.route[i]
            to = evtol.route[i + 1]
            TravelLength = Float32(dist[(from, to)])

            # Charge only the shortfall needed to land at or above bmin,
            # capped by available turnaround time and headroom to bmax.
            need      = BatteryNeeded(TravelLength, battery_per_km)
            shortfall = need - (BatteryLevel[i] - bmin)
            time_poss = BatteryCharged(Float16(evtol.turnaroundTime[i]), ec)
            headroom  = bmax - BatteryLevel[i]
            charge_i  = min(max(0.0f0, shortfall), time_poss, headroom)

            BatteryLevel[i + 1] = BatteryLevel[i] + charge_i - need

            if BatteryLevel[i+1] + need > bmid
                for j in bmid+1:(BatteryLevel[i+1] + need)
                    battery_overrule += b_penalty
                end
            end

            if BatteryLevel[i + 1] < bmin
                feasible = false
            end
        end

        push!(battery_levels, BatteryLevel) 
    end 

    return feasible, battery_levels, battery_overrule
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

function FeasibilityCheck(bmax::Float32, bmid::Float32, bmin::Float32, dist::Dict{Tuple{Int,Int},Float64}, ec::Float32, battery_per_km::Float32,
    evtols::allPlaneSolution,rt::Matrix{Int}, ET::Int, T::Int, V::Int,cap_v::Dict{},b_penalty)

    P = zeros(Int32, 4)

    feasible, battery_levels, battery_overrule = FeasibleBattery(evtols, bmax, bmid, bmin, dist, ec, battery_per_km, b_penalty)

    if feasible == false
        P[1] = 1_000_000
    end

    if feasible == true
        P[1] = battery_overrule 
    end

    if FeasibleCompletionTime(evtols, rt, ET) == false
        P[2] = 1_000_000
    end

    if FeasibleVertiportCapacity(evtols, rt, T, V, cap_v, ET) == false
        P[3] = 1_000_000
    end


    return P
end