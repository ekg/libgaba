#! /usr/bin/env python
# encoding: utf-8

def options(opt):
	opt.load('compiler_c')

	opt.recurse('arch')

def configure(conf):
	conf.load('ar')
	conf.load('compiler_c')

	conf.recurse('arch')

	conf.env.append_value('CFLAGS', '-O3')
	conf.env.append_value('CFLAGS', '-Wall')
	conf.env.append_value('CFLAGS', '-std=c99')
	# conf.env.append_value('CFLAGS', '-fPIC')

	if conf.env.CC_NAME == 'icc':
		# FIXME: dirty hack to pass '-diag-disable remark', current waf does not support space-separated options
		# note: https://groups.google.com/forum/#!topic/waf-users/TJKhe04HGQc
		conf.env.append_value('CFLAGS', '-diag-disable')
		conf.env.append_value('CFLAGS', 'remark')
		conf.env.append_value('CFLAGS', '-inline-forceinline')
	elif conf.env.CC_NAME == 'gcc':
		conf.env.append_value('CFLAGS', '-Wno-unused-variable')
		conf.env.append_value('CFLAGS', '-Wno-unused-but-set-variable')
		conf.env.append_value('CFLAGS', '-Wno-unused-result')
		conf.env.append_value('CFLAGS', '-Wno-format')				# ignore incompatibility between unsigned long long and int64_t
	elif conf.env.CC_NAME == 'clang':
		conf.env.append_value('CFLAGS', '-Wno-unused-variable')
	else:
		pass

	conf.env.append_value('OBJ_GABA', ['gaba.o'])


def build(bld):

	bld.recurse('arch')

	bld.objects(source = 'gaba.c', target = 'gaba.o', includes = ['.'])

	bld.stlib(
		source = ['unittest.c'],
		target = 'gaba',
		includes = '.',
		use = bld.env.OBJ_GABA)

	bld.program(
		source = ['unittest.c'],
		target = 'unittest',
		includes = ['.'],
		use = bld.env.OBJ_GABA,
		defines = ['TEST'])

	bld.program(
		source = ['bench.c'],
		target = 'bench',
		includes = ['.'],
		use = bld.env.OBJ_GABA,
		defines = ['BENCH'])
