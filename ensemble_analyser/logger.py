


def save_snapshot(output, confs):
    with open(output, 'w') as f:
        for i in confs:
            f.write(f'{i}\n')

    return None