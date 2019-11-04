def split(container, count):
	"""
	Simple function splitting a container into equal length chunks.
	Order is not preserved but this is potentially an advantage depending on
	the use case.
	"""
	return [container[_i::count] for _i in range(count)]
