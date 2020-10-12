def substring(string, substringList=[]):
    # Recursive function to create all the substrings
    #   of the given string

    if len(string) == 0:
        return substringList

    else:
        substringList.append(string)
        substring(string[1:], substringList)
        return substringList

print substring("bananas")