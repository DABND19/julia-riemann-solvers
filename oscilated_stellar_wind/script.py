import json


if __name__ == '__main__':
    with open('data/oscilated_stellar_wind_15000.jsonl') as f:
        data = list(map(json.loads, f))

    with open('input.json', 'w') as f:
        json.dump(data[-1], f)
