

T = Float32



X = Vector{Float32}[[-0.053458035, -0.9463467, -2.1219254f-5, -3.4570694f-6, -0.0001705885], [-0.94735324, -0.052405387, -3.5464764f-6, -2.259016f-5, -0.00021523237], [-0.00016996264, -0.00011795759, 3.2782555f-7, 3.2782555f-7, -0.99971277], [-0.0013440549, 0.0013250709, -0.45362252, -0.54635215, -6.377697f-6], [-0.94637513, -0.0535253, 5.066395f-7, -1.8537045f-5, -8.1539154f-5], [-0.042866617, -0.9569653, -0.13819799, 0.13818374, -0.00015375018], [0.0007413626, -0.0011407733, 0.0, 2.9802322f-7, -0.9996009], [-0.052548885, -0.9473598, -2.1517277f-5, -3.4570694f-6, -6.633997f-5], [-0.9472618, -0.05250725, -2.9802322f-8, -1.9550323f-5, -0.00021135807], [-0.052453607, -0.94745517, -1.8388033f-5, -8.34465f-7, -7.200241f-5], [-0.9274003, -0.07248375, -0.1379239, 0.1378937, -8.574128f-5], [-0.9464545, -0.053430766, -3.5464764f-6, -2.2292137f-5, -8.8870525f-5], [-0.95765793, -0.042107612, 0.13788277, -0.1378994, -0.00021776557], [-0.00014439225, -0.00010713935, -2.682209f-7, -2.9802322f-7, -0.9997479], [-0.9472731, -0.052509964, 5.364418f-7, -1.8835068f-5, -0.00019866228], [-0.04216072, -0.9577657, -0.13845155, 0.13843709, -5.9247017f-5]]
#X = collect( randn(T, 5) for _ = 1:20 )

acceptance_factor = convert(T, 0.99) # a value in (0,1). closest to 1 means almost only the part with the maximum dev is added to P.
# A larger value tend to yield fewer parts in the final P.

max_dev = convert(T, 0.2)
level_config = SL.UseMaxDeviation(max_dev)

metric = SL.geteuclideanmetric()
P, iters_ran = SL.iteratedsl(
    level_config,
    SL.geteuclideanmetric(),
    X;
    acceptance_factor = acceptance_factor,
    max_iter = 100, # how many iterations we allow.
)
@show length(P)

P_flat = collect( Iterators.flatten(P))

# P should contain unique integers from 1:length(X)
@assert norm(unique(sort(P_flat)) - collect(1:length(X))) < eps(T)*10


nothing
