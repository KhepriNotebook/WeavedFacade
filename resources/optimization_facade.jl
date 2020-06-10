cross_products_pf(vs::Vecs) =
    if length(vs) == 1
        vxyz(0, 0, 0)
    else
        cross(vs[1], vs[2]) + cross_products_pf(vs[2:end])
    end

normal_polygon_pf(pts) = #usar ::Array{Locs} ?
    let ps=[pts[2:end]..., pts[1]]
        vs=[p2-p1 for (p1,p2) in zip(pts,ps)]
        unitized(cross_products_pf(vs))
    end

transpose_matrix_pf(matrix) = #usar ::Matrix ?
    [[row[i] for row in matrix]
    for i in 1:length(matrix[1])]

normal_square_pf(p0::Loc, p1::Loc, p2::Loc, p3::Loc) =
        normal_polygon_pf([p1 - p0, p2 - p1, p3 - p2, p0 - p3])

normals_surface_pf(ptss)=
    [[normal_square_pf(p0,p1,p2,p3)
     for (p0,p1,p2,p3) in zip(pts[1:end-1],pts[2:end], next_pts[1:end-1], next_pts[2:end])]
    for (pts,next_pts) in zip(ptss[1:end-1],ptss[2:end])]

normals_surface_extra_row_pf(ptss)=
    let ptss1 = [vcat(pts,pts[end]) for pts in normals_surface_pf(ptss)]
      [ptss1..., ptss1[end]]
  end

struct Surf_pf
    f::Function
    u0::Real
    u1::Real
    v0::Real
    v1::Real
end

sub_surf_pf(f::Function, u0::Real, u1::Real, v0::Real, v1::Real)=
    Surf_pf(f, u0, u1, v0, v1)

surf_pts_pf(f::Function, u0::Real, u1::Real, v0::Real, v1::Real, n::Int, m::Int)=
    map_division(f, u0, u1, n, v0, v1, m)

stripe_2D_pts2_pf(pts, amp, bending_pattern, nvs_stripe)=
    [let p1 = pt1 + nv*up_or_down*amp,
         p2 = pt2 + nv*next_up_or_down*amp
      [p1, intermediate_loc(p1, p2, 0.5)]
     end
     for (pt1, pt2, nv, up_or_down, next_up_or_down)
     in zip(pts, [pts[2:end]...,pts[end]], nvs_stripe,
            cycle(bending_pattern),
            cycle([bending_pattern[2:end]...,bending_pattern[1]]))]

stripes_2D_ptss2_pf(ptss, amp, weave_pattern)=
    let nvs = normals_surface_extra_row_pf(ptss)
     [vcat(stripe_2D_pts2_pf(pts, amp, bending_pattern, nv_surf)...)[1:end-1]
     for (pts, nv_surf, bending_pattern) in zip(ptss, nvs, cycle(weave_pattern))]
    end

weave_ptss_centered_pf(s, nu, nv, amp, weave_pattern, stripes_widths, is_verical_straight, is_horizotal_straight)=
    let Δu_ini = maximum(stripes_widths[1])/2,
        Δv_ini = maximum(stripes_widths[2])/2,
        ptss_v = surf_pts_pf(s.f, s.u0 + Δu_ini, s.u1 - Δu_ini, s.v0, s.v1, nu - 1, nv - 1), #corrigir dominios finais!!!
        ptss_u = transpose_matrix_pf(surf_pts_pf(s.f, s.u0, s.u1, s.v0 + Δv_ini, s.v1 - Δv_ini, nu - 1, nv - 1)), #corrigir dominios finais!!!
        inverse_weave_pattern = transpose_matrix_pf(weave_pattern),
        amp_v = is_verical_straight ? 0 : amp,
        amp_u = is_horizotal_straight ? 0 : amp
      [stripes_2D_ptss2_pf(ptss_v, amp_v, weave_pattern),
       stripes_2D_ptss2_pf(ptss_u, amp_u, inverse_weave_pattern)]
   end

# esta sem thicken #
stripe_centered_v_pf(pts, stripe_width, stripe_thick, f_rotations, f_width) =
    thicken(
        surface_grid([map(pt-> pt-vpol(stripe_width/2 + f_width(pt), f_rotations[1](pt)), pts),
                      map(pt-> pt+vpol(stripe_width/2 + f_width(pt), f_rotations[2](pt)), pts)]),
            stripe_thick)

# esta sem thicken #
stripe_centered_u_pf(pts, stripe_width, stripe_thick, f_rotations, f_width) = #permitir a stripe variar a width em certas zonas
    thicken(
        surface_grid([map(pt-> pt-vpol(stripe_width/2 + f_width(pt), f_rotations[1](pt)), pts),
                      map(pt-> pt+vpol(stripe_width/2 + f_width(pt), f_rotations[2](pt)), pts)]),
            stripe_thick)

stripes_with_rotation_pf(ptss, stripes_widths, thickness, f_widths, f_rotations)=
    vcat([width==0 ?
          empty_shape() :
          stripe_centered_v_pf(pt, width, thickness[1], f_rotations[1], f_widths[1])
          for (pt,width) in zip(ptss[1], Iterators.Cycle(stripes_widths[1]))],
         Khepri.with(current_cs, cs_from_o_vx_vy(u0(), vz(1), vy(1))) do
             [width==0 ?
              empty_shape() :
              stripe_centered_u_pf(pt, width, thickness[2], f_rotations[2], f_widths[2])
              for (pt,width) in zip(ptss[2], Iterators.Cycle(stripes_widths[2]))]
         end)

weave_optimize_light_views_pf(s, n_u, amp, weave_pattern, stripes_widths, stripes_thickness,
                           is_vertical_straight, is_horizontal_straight,
                           f_vstripes_width, f_vstripes_rotations)=
    let ptss = weave_ptss_centered_pf(s, n_u, amp, weave_pattern, stripes_widths, is_vertical_straight, is_horizontal_straight)
      stripes_with_rotation_pf(ptss, stripes_widths, stripes_thickness, f_vstripes_width, f_vstripes_rotations)
    end

function facade_pf(; m_stripes_south03=10, angle_rotation=pi/2, f_thicken=pi/2)

    delete_all_shapes()

    f_south_length = 12.55
    f_south_height = 3
    f_length = 10
    amp3 = 0.12
    n_stripes_south03 = ceil(Int, f_length/0.35) # Pera: substituí floor por ceil
    s1= Surf_pf((i,j)-> xyz(i,-0.18,j), -0.3, f_south_length, 0.3, 0.3 + f_south_height)
    stripe_widths_south03 = [vcat(collect(repeated(0.1, 8)),
                                  collect(repeated(0.0, 13)),
                                  collect(repeated(0.1, 8))),
                                  [f_south_height/m_stripes_south03]]
    weave_PLEA3=vcat(collect(take(cycle([[-1,+1,-1,+1],[+1,-1,+1,-1]]), 8)),
                     collect(repeated([-1,-1,-1,-1], 13)),
                     collect(take(cycle([[+1,-1,+1,-1],[-1,+1,-1,+1]]), 8)))
    pts_south03= weave_ptss_centered_pf(s1, n_stripes_south03, m_stripes_south03, amp3, weave_PLEA3, stripe_widths_south03, true, false)

    f_east_length = 9.85
    f_east_height = 3
    n_stripes_east01 = ceil(Int, f_length/0.5) # Pera: substituí floor por ceil
    m_stripes_east01 = m_stripes_south03
    s2= Surf_pf((i,j)-> xyz(f_south_length+0.15,i,j), 0, f_east_length, 0.3, 0.3 + f_east_height)
    stripe_widths_east01 = [[i for i in 0.15:-0.005:0.02], [f_east_height/m_stripes_south03]]
    weave_PLEA_east1=[[-1,+1,-1,+1],
                      [+1,-1,+1,-1],
                      [-1,+1,-1,+1],
                      [+1,-1,+1,-1]]
    pts_east01= weave_ptss_centered_pf(s2, n_stripes_east01, m_stripes_east01, amp3, weave_PLEA_east1, stripe_widths_east01, true, false)
    stripe_widths_east02 = [[i for i in 0.1:0.003:0.2], [f_east_height/m_stripes_east01]]
    pts_east02= weave_ptss_centered_pf(s2, n_stripes_east01, m_stripes_east01, amp3, weave_PLEA_east1, stripe_widths_east02, true, false)

        let ϵ = 0.15
            stripes_with_rotation_pf(pts_south03,
                                  stripe_widths_south03,
                                  [0.01, 0.01],
                                  [i->0, i->0],
                                  [[i->0, i-> 0],
                                   [i->angle_rotation*(ϵ < i.x/f_south_length < 1-ϵ ? sin(pi*(i.x/f_south_length - ϵ)/(1-2ϵ)) : 0),
                                    i->angle_rotation*(ϵ < i.x/f_south_length < 1-ϵ ? sin(pi*(i.x/f_south_length - ϵ)/(1-2ϵ)) : 0)]])

            stripes_with_rotation_pf(pts_east02,
                                  stripe_widths_east02,
                                  [0.01, 0.01],
                                  [i->0, i->-0.15*sin(f_thicken*i.y/f_east_length)*sin(pi*(i.z-0.3)/f_east_height)],
                                  [[i->angle_rotation*(-1+cos(2pi*(i.z-0.3)/f_east_height))/2.5, i->angle_rotation*(-1+cos(2pi*(i.z-0.3)/f_east_height))/2.5],
                                   [i->0, i->0]])
                           end
end
