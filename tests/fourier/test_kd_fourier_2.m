clear
clf

[pos,els]=load_pdetool_geom( 'square', 'numrefine', 3 );

k1 = 3; 
k2 = 7;

w_k = [
    k1, k2
    k1, k2
    k1, k2
    k1, k2
    ]/2;

p_k = [
    0, 0
    1, 1
    0, 1
    1, 0
    ] * 0.25;

TB = { w_k, p_k };

for a=linspace(0,1,41)
%     u1 = trig_eval([1 -1 0  0], TB, pos);
%     u0 = trig_eval([0  0 1 -1], TB, pos);
%     u0 = trig_eval([1 1 0  0], TB, pos);
%     u0 = trig_eval([0  0 1 1], TB, pos);

    u0 = trig_eval([0 -1 0  0], TB, pos);
    u1 = trig_eval([1  0 0  0], TB, pos);
    u = (1-a)*u0 + a * u1;
    plot_field(pos, els, u', 'show_mesh', false)
    drawnow
end
