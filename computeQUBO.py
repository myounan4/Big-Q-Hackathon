import numpy as np
import matplotlib.pyplot as plt
from pulser import Pulse, Sequence, Register
from pulser_simulation import QutipEmulator
from pulser.devices import Chadoq2
from pulser.waveforms import InterpolatedWaveform
from scipy.optimize import minimize
from scipy.spatial.distance import pdist, squareform

def top_n_bitstrings(count_dict, n):
    # Find the top n bit strings with the highest counts
    top_strings = [item[0] for item in count_dict.most_common(n)]
    
    # Convert each bit string to a list of integers
    return [[int(bit) for bit in string] for string in top_strings]

def computeQUBO(Q, delta_params, amp_params):
    bitstrings = [np.binary_repr(i, len(Q)) for i in range(2 ** len(Q))]
    costs = []
    # this takes exponential time with the dimension of the QUBO
    for b in bitstrings:
        z = np.array(list(b), dtype=int)
        cost = z.T @ Q @ z
        costs.append(cost)
    zipped = zip(bitstrings, costs)
    sort_zipped = sorted(zipped, key=lambda x: x[1])
    print(sort_zipped[:3])
    def evaluate_mapping(new_coords, *args):
        """Cost function to minimize. Ideally, the pairwise
        distances are conserved"""
        Q, shape = args
        new_coords = np.reshape(new_coords, shape)
        new_Q = squareform(Chadoq2.interaction_coeff / pdist(new_coords) ** 6)
        return np.linalg.norm(new_Q - Q)
    shape = (len(Q), 2)
    costs = []
    np.random.seed(0)
    x0 = np.random.random(shape).flatten()
    res = minimize(
        evaluate_mapping,
        x0,
        args=(Q, shape),
        method="Nelder-Mead",
        tol=1e-6,
        options={"maxiter": 200000, "maxfev": None},
    )
    coords = np.reshape(res.x, (len(Q), 2))
    qubits = dict(enumerate(coords))
    reg = Register(qubits)

    #reg.draw(
    #    blockade_radius=Chadoq2.rydberg_blockade_radius(1.0),
    #    draw_graph=False,
    #    draw_half_radius=True,
    #)

    LAYERS = 2
    
    # Parametrized sequence
    seq = Sequence(reg, Chadoq2)
    seq.declare_channel("ch0", "rydberg_global")

    t_list = seq.declare_variable("t_list", size=LAYERS)
    s_list = seq.declare_variable("s_list", size=LAYERS)

    for t, s in zip(t_list, s_list):
        pulse_1 = Pulse.ConstantPulse(1000 * t, 1.0, 0.0, 0)
        pulse_2 = Pulse.ConstantPulse(1000 * s, 0.0, 1.0, 0)

        seq.add(pulse_1, "ch0")
        seq.add(pulse_2, "ch0")

    seq.measure("ground-rydberg")
    def quantum_loop(parameters):
        params = np.array(parameters)
        t_params, s_params = np.reshape(params.astype(int), (2, LAYERS))
        assigned_seq = seq.build(t_list=t_params, s_list=s_params)
        simul = QutipEmulator.from_sequence(assigned_seq, sampling_rate=0.01)
        results = simul.run()
        count_dict = results.sample_final_state()  # sample from the state vector
        return count_dict
    np.random.seed(123)  # ensures reproducibility of the tutorial
    guess = {
        "t": np.random.uniform(8, 10, LAYERS),
        "s": np.random.uniform(1, 3, LAYERS),
    }
    example_dict = quantum_loop(np.r_[guess["t"], guess["s"]])

    def plot_distribution(C):
        C = dict(sorted(C.items(), key=lambda item: item[1], reverse=True))
        indexes = ["01011", "00111"]  # QUBO solutions
        color_dict = {key: "r" if key in indexes else "g" for key in C}
        plt.figure(figsize=(12, 6))
        plt.xlabel("bitstrings")
        plt.ylabel("counts")
        plt.bar(C.keys(), C.values(), width=0.5, color=color_dict.values())
        plt.xticks(rotation="vertical")
        plt.show()
    #plot_distribution(example_dict)
    
    def get_cost_colouring(bitstring, Q):
        z = np.array(list(bitstring), dtype=int)
        cost = z.T @ Q @ z
        return cost

    def get_cost(counter, Q):
        cost = sum(counter[key] * get_cost_colouring(key, Q) for key in counter)
        return cost / sum(counter.values())  # Divide by total samples
    def func(param, *args):
        Q = args[0]
        C = quantum_loop(param)
        cost = get_cost(C, Q)
        return cost
    scores = []
    params = []

    for repetition in range(20):
        guess = {
            "t": np.random.uniform(1, 10, LAYERS),
            "s": np.random.uniform(1, 10, LAYERS),
        }

        try:
            res = minimize(
                func,
                args=Q,
                x0=np.r_[guess["t"], guess["s"]],
                method="Nelder-Mead",
                tol=1e-5,
                options={"maxiter": 10},
            )
            scores.append(res.fun)
            params.append(res.x)
        except Exception as e:
            pass

    optimal_count_dict = quantum_loop(params[np.argmin(scores)])
    plot_distribution(optimal_count_dict) 
    # We choose a median value between the min and the max
    Omega = np.median(Q[Q > 0].flatten())
    delta_0 = -5  # just has to be negative
    delta_f = -delta_0  # just has to be positive
    T = 4000  # time in ns, we choose a time long enough to ensure the propagation of information in the system
    adiabatic_pulse = Pulse(
        InterpolatedWaveform(T, [1e-9, Omega, 1e-9]),
        InterpolatedWaveform(T, [delta_0, 0, delta_f]),
        0,
    )

    seq = Sequence(reg, Chadoq2)
    seq.declare_channel("ising", "rydberg_global")
    seq.add(adiabatic_pulse, "ising")
    #seq.draw()
    simul = QutipEmulator.from_sequence(seq)
    results = simul.run()
    final = results.get_final_state()
    count_dict = results.sample_final_state()
    plot_distribution(count_dict)
    
    
    cost = []

    for T in 1000 * np.linspace(1, 10, 10):
        seq = Sequence(reg, Chadoq2)
        seq.declare_channel("ising", "rydberg_global")
        adiabatic_pulse = Pulse(
            InterpolatedWaveform(T, [1e-9, Omega, 1e-9]),
            InterpolatedWaveform(T, [delta_0, 0, delta_f]),
            0,
        )
        seq.add(adiabatic_pulse, "ising")
        simul = QutipEmulator.from_sequence(seq)
        results = simul.run()
        final = results.get_final_state()
        count_dict = results.sample_final_state()
        cost.append(get_cost(count_dict, Q) / 3)
        plt.figure(figsize=(12, 6))
    plt.plot(range(1, 11), np.array(cost), "--o")
    plt.xlabel("total time evolution (µs)", fontsize=14)
    plt.ylabel("cost", fontsize=14)
    plt.show()
    
    n = 10 # number of bit strings we want 
    
    return top_n_bitstrings(count_dict, n)
