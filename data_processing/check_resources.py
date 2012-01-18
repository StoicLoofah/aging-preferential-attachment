f = open('PINTS03-06', 'r')
resources = {}
lines_read = 0
for line in f:
    content = line.split('\t')
    user = int(content[1])
    resource = int(content[2])
    if resource in resources:
        if user not in resources[resource]:
            print('FOUND')
            print(resource)
            print(resources[resource])
            print(user)
            break
        else:
            resources[resource].add(user)
    else:
        resources[resource] = set([user])
    lines_read += 1
    if lines_read % 1000 == 0:
        print(lines_read)
f.close()
