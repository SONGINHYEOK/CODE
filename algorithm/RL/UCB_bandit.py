import numpy as np

class UCB_Bandit:
    def __init__(self, arms=10):
        self.arms = arms
        self.q_values = np.zeros(arms)  # Estimated values of each arm
        self.action_counts = np.zeros(arms)  # Number of times each arm has been pulled
        self.total_pulls = 0

    def ucb(self, arm):
        if self.action_counts[arm] == 0:
            return float('inf')  # Return infinity for unexplored arms
        else:
            exploration_term = np.sqrt(2 * np.log(self.total_pulls) / self.action_counts[arm])
            return self.q_values[arm] + exploration_term

    def select_arm(self):
        ucb_values = [self.ucb(arm) for arm in range(self.arms)]
        return np.argmax(ucb_values)

    def play(self, rounds=1000):
        for _ in range(rounds):
            arm = self.select_arm()
            reward = self.play_arm(arm)
            self.update_values(arm, reward)

    def play_arm(self, arm):
        # Simulate pulling the selected arm and getting a reward (1 or 0, for example)
        # Here, you can replace this with your actual bandit play logic
        rate = np.random.rand()  # Simulating a random success rate for the selected arm
        reward = 1 if rate > 0.5 else 0
        return reward

    def update_values(self, arm, reward):
        self.action_counts[arm] += 1
        self.total_pulls += 1
        # Update the estimated value of the selected arm using a simple average
        self.q_values[arm] += (reward - self.q_values[arm]) / self.action_counts[arm]

# Example usage
ucb_bandit = UCB_Bandit(arms=10)
ucb_bandit.play(rounds=1000)

# Print the estimated values for each arm after playing
print("Estimated values:", ucb_bandit.q_values)
