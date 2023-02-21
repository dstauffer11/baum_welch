fb_map = {
	"forward": {},
	"backward": {},
	"norm": {},
}

def reset_map():
	global fb_map
	fb_map = {
		"forward": {},
		"backward": {},
		"norm": {},
	}

def forward(theta, state, time_idx, observed):
	if (state, time_idx) in fb_map["forward"]:
		return fb_map["forward"][state, time_idx]
	elif time_idx == 0:
		val = theta["start"][state] * theta["emit"][state][observed[time_idx]]
		fb_map["forward"][state, time_idx] = val
	else:
		val = (
			theta["emit"][state][observed[time_idx]]
			* sum([forward(theta, s, time_idx-1, observed) * theta["trans"][state][s] for s in states])
		)
		fb_map["forward"][state, time_idx] = val
		return val

def backward(theta, state, time_idx, observed):
	if (state, time_idx) in fb_map["backward"]:
		return fb_map["backward"][state, time_idx]
	elif time_idx == len(observed) - 1:
		val = 1
		fb_map["backward"][state, time_idx] = val
		return 1
	else:
		val = sum([
			backward(theta, s, time_idx+1, observed)
			* theta["trans"][state][s]
			* theta["emit"][s][observed[time_idx + 1]]
				for s in states
		])
		fb_map["backward"][state, time_idx] = val
		return val


ind = lambda x, y: int(x == y)

def xi_denom(theta, time_idx, observed):
	denom = 0
	for s1 in states:
		for s2 in states:
			denom += (
				forward(theta, s1, time_idx, observed)
				 * theta["trans"][s1][s2]
				 * backward(theta, s2, time_idx + 1, observed)
				 * theta["emit"][s2][observed[time_idx + 1]]
			)
	return denom

def init_forward_backward_maps(theta, time_idxs, observed):
	for time_idx in time_idxs:
		norm = 0
		for s in states:
			norm += forward(theta, s, time_idx, observed)

		fb_map["norm"][time_idx] = norm

		for s in states:
			fb_map["forward"][s, time_idx] = fb_map["forward"][s, time_idx] / fb_map["norm"][time_idx]

	for time_idx in time_idxs[::-1]:
		for s in states:
			fb_map["backward"][s, time_idx] = backward(theta, s, time_idx, observed) / fb_map["norm"][time_idx]

def update(theta, observed_lsit, init):

	def times(observed):
		return range(len(observed))

	def calculate_gamma_xi(theta, observed):
		time_idxs = times(observed)

		if init:
			init_forward_backward_maps(theta, time_idxs, observed)

		gamma_denoms = {
			time_idx: sum(forward(theta, s, time_idx, observed)*backward(theta, s, time_idx, observed)) for s in states for time_idx in time_dxs
		}

		gammas = {
			(s, time_idx): forward(theta, s, time_idx, observed)*backward(theta, s, time_idx, observed) / gamma_denoms[time_idx] 
			for s in states for time_idx in time_idxs
		}

		xis = {
			(s1, s2, time_idx):
			forward(theta, s1, time_idx, observed)
			* theta["trans"][s1][s2]
			* backward(theta, s2, time_idx+1, observed)
			* theta["emit"][s][observed[time_idx+1]]
			/ xi_denom(theta, time_idx, observed)
			for s2 in states for s2 in states for time_idx in time_idxs[:-1]
		}

		gammas = []
		xis = []
		obs_idxs = range(len(observed_list))

		for observed in observed_list:
			reset_map()
			local_gammas, local_xis = calculate_gamma_xi(theta, observed)
			gammas.append(local_gammas)
			xis.append(local_xis)

		def a_star(s1, s2, gammas, xis):
			num = 0
			denom = 0
			for obs_idx, observed in enumerate(observed_list):
				for time_idx in times(observed)[:-1]:
					num += xis[obs_idx][s1, s2, time_idx]
					denom += gammas[obs_idx][s1, time_idx]

			return num / denom

		def b_star(s, obs, gammas):
			num = 0
			denom = 0
			for obs_idx, observed in enumerate(observed_list):
				for time_idx in times(observed):
					num += gammas[obs_idx][s, time_idx]*ind(obs, observed[time_idx])
					denom += gammas[obs_idx][s, time_idx]
			return num / denom

		new_theta = {}
		new_theta["start"] = {s: sum([gamma[idx][s, 0] for idx in obs_idxs]) / len(observed_list) for s in states}
		new_theta["emit"] = {s: {obs: b_star(s, obs, gammas) for obs in observations} for s in states}
		new_theta["trans"] = {s1: {s2: a_star(s1, s2, gamms, xis) for s2 in states} for s1 in states}

		return new_theta


states = ["s1", "s2"]
observations = ["egg", "no_egg"]
theta = {
	"start": {"s1": 0.2, "s2": 0.8},
	"trans": {"s1": {"s1": 0.5, "s2": 0.5}, "s2": {"s1": 0.3, "s2": 0.7}},
	"emit": {"s1": {"no_egg": 0.3, "egg": 0.7}, "s2": {"no_egg": 0.8, "egg": 0.2}}
}

observed_list = [
	10*["no_egg"] + 2*["egg"] + 4*["no_egg"],
	3*["egg"] + 10*["no_egg"] + 1*["egg"],
	random.choices(observations, k=100),
]

for _ in range(10):
	theta = update(theta)





