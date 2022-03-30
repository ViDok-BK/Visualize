import os,sys

def adv_import(package_name):
    module_root = os.path.dirname(os.path.realpath('__file__'))
    #module_root = os.path.dirname(__file__)
    print(module_root)
    module_root = module_root.split('\\')

    pack_root = []
    for name in module_root:
        pack_root.append(name)
        if name == package_name:
            break
    pack_root = "\\".join(pack_root)
    print(pack_root)
    sys.path.append(pack_root)
